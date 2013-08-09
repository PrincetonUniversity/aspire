function [G, Z, W, y, out] = cryoEMADM(NN, G, lambda, S, opts)

%------------------------------------------------------------------
% ADMM for Cryo-EM SDP
%
% ctype = 1:
% min  -<S, G> + lambda ||G||_2
% s.t. A(G) = b, G psd
%
% ctype = 0:
% min  -<S, G>
% s.t. A(G) = b, G psd
%      ||G||_2 <= lambda
%
% Author: Zaiwen Wen
% date: 7/14, 2012
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
if ~isfield(opts, 'max_rankZ'); opts.max_rankZ = max(6, floor(NN/4)); end
if ~isfield(opts, 'max_rankW'); opts.max_rankW = max(6, floor(NN/4)); end
if ~isfield(opts, 'eig_ratio_leading'); opts.eig_ratio_leading = 0.01; end

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
nrmb = norm(b);
nrmS = norm(S,inf);

% set up linear constraints
idx = 1:n; idx2 = 1:NN;
col = [NN:n-1]*n+idx2; col = col';

% intial y
AS = ComputeAX(S);
resi = ComputeAX(G) - b;  
W = eye(n);
Z = W;

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
    
    % compute Z
    ATy = ComputeATy(y);
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
        [U, pi] = eigs(B,kk,'lm',opteigs); pi = diag(pi);
    end
    
    if ctype == 1
        zz = projsplx(abs(pi)/lambda);
        zz = (lambda*sign(pi)).*zz;
    elseif ctype == 0
        zz = sign(pi).*max(abs(pi)-lambda/mu, 0);
    end
    
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
    if ctype == 1
        spG = eigs(G, 1, 'lm');
        pobj = -full(sum(sum(S.*G))) + lambda*spG;
        dobj = full(b'*y);
    elseif ctype == 0
        pobj = -full(sum(sum(S.*G)));
        dobj = full(b'*y)-lambda*sum(abs(zz));
    end
        
    gap  = abs(dobj-pobj)/(1+abs(dobj)+abs(pobj));
    resi = ComputeAX(G) - b;
    pinf = norm(resi)/nrmb;
    dinf = norm(S+W+ATy-Z,'fro')/nrmS;

    ainf = max([pinf,dinf,gap]);
    dtmp = pinf/dinf;

    if record >= 1
        %fprintf('%4d   %3.2e   %+3.2e   %+3.2e   %+3.2e   %3.2e   %+3.2e  %4d (%d) (%3.2e,  %3.2e)   %d  %d %d\n', ...
        %    itr, mu,  pobj,  dobj, gap, pinf, dinf,  nev, TcompEigS, ref_inf, ainf, itr_stag,  itmu_pinf, itmu_dinf)%pinf/dinf, atol);

        fprintf('%4d   %3.2e   %+3.2e   %+3.2e   %+3.2e   %3.2e   %+3.2e   %+3.2e  %d  %d  %d %d\n', ...
                itr, mu,  pobj,  dobj, gap, pinf, dinf,  dtmp, kk, nev, itmu_pinf, itmu_dinf);
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



out.pobj = pobj;
out.dobj = dobj;
out.gap  = gap;
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
end

