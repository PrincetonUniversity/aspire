function [G, Z, W, y, out] = cryoEMSDPL12N(NN, G, lambda, C, opts)

%------------------------------------------------------------------
% ADMM for Cryo-EM SDP
%
% min  sum_{i<j} ||c_ij - G_ij c_ji||
% s.t. A(G) = b, G psd
%      ||G||_2 <= lambda
%
% Author: Zaiwen Wen, Lanhui Wang
% date: 7/18, 2012
%------------------------------------------------------------------

%------------------------------------------------------------------
% set up parameters
if ~isfield(opts, 'tol');    opts.tol = 1e-3;   end
if ~isfield(opts, 'mu');     opts.mu = 1;       end
if ~isfield(opts, 'gam');    opts.gam = 1.618;  end
if ~isfield(opts, 'EPS');    opts.EPS = 1e-12;  end
if ~isfield(opts, 'maxit');  opts.maxit = 1000; end
if ~isfield(opts, 'record'); opts.record = 0;   end

% ctype = 1, ||G||_2 in the obj. ctype = 0, in the constraints
if ~isfield(opts, 'ctype');  opts.ctype = 0;    end


if ~isfield(opts, 'adp_proj');  opts.adp_proj = 1; end %1 or 0
if ~isfield(opts, 'max_rankZ'); opts.max_rankZ = max(6, floor(NN/2)); end
if ~isfield(opts, 'max_rankW'); opts.max_rankW = max(6, floor(NN/2)); end

tol     = opts.tol;
mu      = opts.mu;
gam     = opts.gam;
EPS     = opts.EPS;
maxit   = opts.maxit;
record  = opts.record;
ctype   = opts.ctype;

adp_proj = opts.adp_proj;
max_rankZ = opts.max_rankZ;
max_rankW = opts.max_rankW;


% parameters for adjusting mu
if ~isfield(opts, 'adp_mu');        opts.adp_mu = 1;        end %1 or 0
if ~isfield(opts, 'dec_mu');        opts.dec_mu = 0.5;      end
if ~isfield(opts, 'inc_mu');        opts.inc_mu = 2;        end
if ~isfield(opts, 'mu_min');        opts.mu_min = 1e-4;     end
if ~isfield(opts, 'mu_max');        opts.mu_max = 1e4;      end
if ~isfield(opts, 'min_mu_itr');    opts.min_mu_itr = 5;    end
if ~isfield(opts, 'max_mu_itr');    opts.max_mu_itr = 20;   end
if ~isfield(opts, 'delta_mu_l');    opts.delta_mu_l = 0.1;  end
if ~isfield(opts, 'delta_mu_u');    opts.delta_mu_u = 10;   end

adp_mu  = opts.adp_mu;
dec_mu  = opts.dec_mu;  
inc_mu  = opts.inc_mu;  
mu_min  = opts.mu_min;   
mu_max  = opts.mu_max;
max_mu_itr = opts.max_mu_itr;
delta_mu_l = opts.delta_mu_l;
delta_mu_u = opts.delta_mu_u;

itmu_pinf = 0;  
itmu_dinf = 0;

%------------------------------------------------------------------
% set up the SDP problem
n = 2*NN;
m = 3*NN;
b = [ones(n,1); zeros(NN,1)];

% set up linear constraints
idx = 1:n; 
col = ([2:2:n]'-1)*n+[1:2:n]';

% intial values
W = eye(n);
Z = W;

% provide an initial theta here:
Phi = G/mu;
% theta=C2theta(Phi,C,mu);
% % then compute S = Q(theta)
% S = Qtheta(theta,C);
[S, theta]=Qtheta(Phi,C,mu);
S=(S+S')/2;
AS = ComputeAX(S);

resi = ComputeAX(G) - b;

kk = 0; nev = 0;
opteigs.maxit = 100; opteigs.issym = 1;
opteigs.isreal = 1; opteigs.disp = 0; opteigs.tol = 1e-6;

%------------------------------------------------------------------
if record >= 1
    fprintf('%4s %8s %10s %10s %10s %10s %10s\n', ...
        'itr', 'mu', 'pobj', 'dobj', 'gap', 'pinf', 'dinf');
end

% main routine
for itr = 1:maxit
    %------------------------------------------------------------------
    % compute y
    y = -(AS + ComputeAX(W) - ComputeAX(Z)) - resi/mu;
    
    %******************************
    % compute theta
    ATy = ComputeATy(y);
    Phi = W + ATy -Z + G/mu;
    %     theta = C2theta(Phi,C,mu);
    %     % then compute S = Q(theta)
    %     S = Qtheta(theta,C);
    [S, theta]=Qtheta(Phi,C,mu);
    S=(S+S')/2;
    %******************************
    
    % compute Z
    B = S + W + ATy + G/mu; B = (B + B')/2;

    if adp_proj == 0
        %[U, pi] = mexeig(B); %pi = diag(pi);
        [U, pi] = eig(B); pi = diag(pi);
    elseif adp_proj == 1
        if itr == 1
            kk = max_rankZ;
        else
            if  kk > 0
                drops = zz(1:end-1)./zz(2:end);
                [dmx,imx] = max(drops);
                rel_drp = (nev-1)*dmx/(sum(drops)-dmx);
                %rel_drp
                if rel_drp > 10
                    kk = max(imx, 6);
                else
                    kk = kk + 3;
                end
            else
                kk = 6;
            end
        end
        kk = min(kk, n);
        [U, pi] = eigs(B,kk,'lm',opteigs); pi = diag(pi);
    end
    
    zz = sign(pi).*max(abs(pi)-lambda/mu, 0);
    nD = abs(zz) > 0;   kk = nnz(nD);  
    if kk > 0;
        zz = zz(nD);
        Z = U(:,nD)*diag(zz)*U(:,nD)';
    else
        Z = zeros(n);
    end
    
    % compute W
    H = Z - S - ATy - G/mu;     H = (H + H')/2;
    
    if adp_proj == 0
        [V,D] = eig(H);
        D = diag(D);
        nD = D>EPS;
        nev = nnz(nD);
        if nev < n/2 % few positive eigenvalues
            W = V(:,nD)*diag(D(nD))*V(:,nD)';
            WmH = W - H;
        else  % few negative eigenvalues
            nD = ~nD;
            WmH = V(:,nD)*diag(-D(nD))*V(:,nD)';
            W = WmH + H;
        end
    elseif adp_proj == 1 
        % compute G = W - H since G is expected to be rank 3
        % estimate rank, 6 is a safeguard here
        if itr == 1
            nev = max_rankW;
        else
            if  nev > 0
                %nev
                drops = dH(1:end-1)./dH(2:end);
                [dmx,imx] = max(drops);  
                rel_drp = (nev-1)*dmx/(sum(drops)-dmx);
                %rel_drp
                if rel_drp > 50
                    nev = max(imx, 6); 
                else
                    nev = nev + 5;
                end
                %nev
            else
                nev = 6;
            end
        end
        
        % computation of W and H
        nev = min(nev, n);
        [V,dH] = eigs(-H, nev,'la', opteigs);
        dH = diag(dH); nD = dH>EPS; nev = nnz(nD); 
        %dH'
        if nev > 0
            dH = dH(nD);
            WmH = V(:,nD)*diag(dH)*V(:,nD)';
            W = WmH + H;
        else
            WmH = sparse(n,n);
            W = H;
        end
    end
    
    % update G
    %G = (1-gam)*G + gam*mu*(W-H);
    G = (1-gam)*G + gam*mu*WmH;
    
    %------------------------------------------------------------------
    % check optimality
    %pobj = -full(sum(sum(S.*G)));
    %dobj = full(b'*y)-lambda*sum(abszz);
    %gap  = abs(dobj-pobj)/(1+abs(dobj)+abs(pobj));
   
    resi = ComputeAX(G) - b;
    
    nrmb = max(norm(b),1);
    nrmS = max(norm(S,inf),1);
    
    pinf = norm(resi)/nrmb;
    dinf = norm(S+W+ATy-Z,'fro')/nrmS;

    %ainf = max([pinf,dinf,gap]);
    ainf = max([pinf,dinf]);
    dtmp = pinf/dinf;

    if record >= 1
        %fprintf('%4d   %3.2e   %+3.2e   %+3.2e   %+3.2e   %3.2e   %+3.2e   %+3.2e  %d  %d  %d %d\n', ...
        %        itr, mu,  pobj,  dobj, gap, pinf, dinf,  dtmp, kk, nev, itmu_pinf, itmu_dinf);
    
        fprintf('%4d   %3.2e   %3.2e   %+3.2e   %+3.2e  %d  %d  %d %d\n', ...
                itr, mu,  pinf, dinf,  dtmp, kk, nev, itmu_pinf, itmu_dinf);
    
    end

    if ainf <= tol
        out.exit = 'optimal';
        break;
    end
    

    % update mu adpatively
    if adp_mu == 1
        if (dtmp <= delta_mu_l)
            itmu_pinf = itmu_pinf + 1;  itmu_dinf = 0;
            if itmu_pinf > max_mu_itr
                %mu = max(mu*dec_mu, mu_min); itmu_pinf = 0;
                mu = max(mu*inc_mu, mu_min); itmu_pinf = 0;
            end
        elseif dtmp > delta_mu_u
            itmu_dinf = itmu_dinf + 1;  itmu_pinf = 0;
            if itmu_dinf > max_mu_itr
                %mu = min(mu*inc_mu, mu_max); itmu_dinf = 0;
                mu = min(mu*dec_mu, mu_max); itmu_dinf = 0;
            end
        end
    end
  
end



%out.pobj = pobj;
% out.dobj = dobj;
% out.gap  = gap;
out.itr = itr;
out.pinf = pinf;
out.dinf = dinf;
% End main routine
%--------------------------------------------------------------------------


	% compute A'*y 
    function ATy = ComputeATy(y)
        ATy = sparse(n,n);
        ATy(col) = (sqrt(2)/2)*y(n+1:m);
        ATy = ATy + ATy'+ sparse(idx, idx, y(1:n), n, n);
    end

    % compute A*X
    function AX = ComputeAX(X)
        AX = [spdiags(X,0); sqrt(2)*X(col)];
    end

%     function theta = C2theta(Phi, C, mu)
%         % Update theta
%         %
%         % Lanhui Wang
%         % Jul 18, 2012
%         K = size(C,2);
%         theta = zeros(2,K,K);
%         
%         for k1 = 1:K
%             for k2 = k1+1:K
%                 c_ij = C(:,k1,k2);
%                 c_ji = C(:,k2,k1);
%                 phi = Phi(2*k1-1:2*k1,2*k2-1:2*k2);
%                 theta(:,k1,k2) = (c_ij-mu*phi*c_ji)/norm(mu*phi*c_ji-c_ij);
%             end
%         end
%         
%     end
% 
% 
%     function S = Qtheta(theta,C)
%         % Update S
%         %
%         % Lanhui Wang
%         % Jul 18, 2012
%         
%         K=size(C,2);
%         S=zeros(2*K);
%         
%         for k1 = 1:K
%             for k2 = k1+1:K
%                 theta_ij = theta(:,k1,k2);
%                 c_ji = C(:,k2,k1);
%                 S(2*k1-1:2*k1,2*k2-1:2*k2) = theta_ij*c_ji';
%             end
%         end
%         
%         S=(S+S')/2;
%     end



end

