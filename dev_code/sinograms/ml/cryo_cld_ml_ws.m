function [S,Z,Zi,PCs,mse,shifts] = cryo_cld_ml_ws(X,L,K,KNN,max_itr,noise_psd,verbose,max_shift,shift_step,Xclean)%% cryo_cld_ml Common lines detection in cryo-em images using EM.
%% cryo_cld_ml Common lines detection in cryo-em images using EM.
%
% Syntax:  [S,Z,Zi,PCs] = cld_ml(X,L,K,max_itr,noise_psd,verbose)
%
% Inputs:
%    X - matrix of noisy polar fourier lines (vectors).
%    L - angular resolution.
%    K - number of clusters.
%    KNN - number of clusters for the Estep (KNN<=K).
%    max_itr - maximum number of iterations (maybe to be substitute with 
%              a stoping condition...).
%    noise_psd - noise power spectrum matrix.
%    max_shift - maximum shift (in pixels).
%    shift_step - shift step (in pixels).
%
% Outputs:
%    S - clusters centers.
%    Z - likelihoods.
%    Zi - cluster indicators.
%    PCs - principal componenets of X.
%    mse - MSE in S (for debug)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Mor Cohen 


narginchk(3, 10);
mse_check=0;
switch nargin
    case 3
        KNN=50;
        max_itr=0;
        noise_psd=1;
        verbose=1;      
    case 4
        max_itr=0;
        noise_psd=1;
        verbose=1;   
    case 5
        noise_psd=1;
        verbose=1;   
    case 6
        verbose=1;   
    case 8
        mse_check=1;        
        if verbose
            fprintf('debug mode\n');
        end
end

if ~isinteger(K)
    K=ceil(K);
end

N=size(X,2)/L;
if K<N
    error('K input must be greater or equal N.')
end
%NL=N*L;
rmax=size(X,1);
rk2=(0:(rmax-1))';
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end

C = pinv(noise_psd,0.01);
%% intialization (using PCA)
Z=0;
Zi=0;
[S,h,PCs,init_pick]=init(X,N,L,K);

S_stack=zeros(size(S,1),size(S,2)*n_shifts);


mse=zeros(1,max_itr+2);
if mse_check
    Sclean=Xclean(:,init_pick);
    Sclean=cryo_raynormalize(Sclean);
    Y=cryo_raynormalize(X(:,init_pick));
    d=abs(Y-Sclean).^2;
    mse(1)=mean(d(:));
    d=abs(S-Sclean).^2;
    mse(2)=mean(d(:));
    if verbose
        msg=sprintf(' mse=%f\n',mse(2));
        fprintf('%s',msg);
    end
end
%% main loop
for itr=1:max_itr
    
    t1=clock;
    
    if verbose
        msg=sprintf('itr=%2d/%2d ',itr,max_itr);
        fprintf('%s',msg);
    end
    
    if verbose
        fprintf('E step');
    end
    % Instead of processing the shifts in a loop, one shift at a time, stack
    % all shifts of a into each cluster
    for k=1:n_shifts
        S_stack(:,(k-1)*K+1:k*K)=bsxfun(@times,S,shift_phases(:,k));
    end
    
    [Z,Zi]=Estep(sqrt(C)*X,L,h,sqrt(C)*S_stack,KNN);
    % map stack index to cluster and shift index
    [Zi,shifts]=ind2sub([K n_shifts],Zi);
    if verbose
        fprintf('/M step ');
    end 
    
    S=Mstep_S(X,Z,Zi,S,h,K,KNN,shifts,max_shift,shift_step);
    h=Mstep_h(X,Z,Zi,S,h,C,N,L);
    
    t2=clock;
    t=etime(t2,t1);
    
    if verbose
        msg=sprintf(' t=%7.5f\n',t);
        fprintf('%s',msg);
    end
    
    if mse_check
        S=cryo_raynormalize(S);
        d=abs(S-Sclean).^2;
        mse(2+itr)=mean(d(:));
        if verbose
            msg=sprintf(' mse=%f\n',mse(2+itr));
            fprintf('%s',msg);
        end
    end
end
%if max_itr>0
%    Z=Z(1,:);
%    Zi=Zi(1,:);
%end

function [Z,Zi]=Estep(X,L,ht,St,KNN)
NL=size(X,2);
N=ceil(NL/L);
Zi=zeros(KNN,L,N); % sliced variable
Z=zeros(KNN,L,N); % sliced variable
parfor n=1:N
    nl=(n-1)*L+(1:L);
    [Zi(:,:,n),Z(:,:,n)]=knn_opt(bsxfun(@times,X(:,nl),1./ht(nl)'),St,KNN);
end

Zi=reshape(Zi,KNN,NL);
Z=reshape(Z,KNN,NL);
Z=bsxfun(@times,Z,(ht').^2); % multiply back by h^2 
Z=exp(-Z);

function [h]=Mstep_h(X,Z,Zi,St,h,C,N,L)

% tic;
% for n=1:N*L
%     Sk=St(:,Zi(:,n));
%     Znlk=Z(:,n);
%     Xnl=X(:,n);
%     % XXX Why the real is only on the nominator? Why real at all?
%     h(n)=real((C*Xnl)'*Sk)*Znlk/(Znlk'*diag(Sk'*(C*Sk)));
% end
% toc

CSt=C*St;
dn1=sum(conj(St).*CSt); dn1=dn1.';
CX=C*X; CX=CX';
parfor n=1:N*L
    Znlk=Z(:,n);
    denom=sum(conj(Znlk).*dn1(Zi(:,n),:));
    h(n)=real(CX(n,:)*St(:,Zi(:,n)))*Znlk/denom;
end

%assert(norm(h(:)-h2(:))/norm(h(:))<1.0e-13);


function S=Mstep_S(X,Z,Zi,S,ht,K,KNN,shifts,max_shift,shift_step) 

rmax=size(X,1);
rk2=(0:(rmax-1))';
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end

% for k=1:K
%     ind=find(Zi==k);
%     if isempty(ind)
%         continue;
%     end
%     nl=ceil(ind/KNN);
%     Unlk=shift_phases(:,shifts(ind));  
%     Znlk=Z(ind);
%     S(:,k)=sum(bsxfun(@times,(ht(nl).*Znlk)',Unlk.*X(:,nl)),2)/sum(ht(nl).^2.*Znlk);
% end
[Zis,ZiI]=sort(Zi(:));
% tic;
% for k=1:K
%     [lower_index,upper_index] = bsearch(Zis,k,k);
%     if isempty(lower_index) || isempty(upper_index)
%         continue;
%     end
%     ind=ZiI(lower_index:upper_index);
%         nl=ceil(ind/KNN);
%     Unlk=shift_phases(:,shifts(ind));  
%     Znlk=Z(ind);
%     
%     % XXX Why is ht of length N*L and not of length N?
%     % XXX It seems that you access
%     S(:,k)=sum(bsxfun(@times,(ht(nl).*Znlk)',Unlk.*X(:,nl)),2)/sum(ht(nl).^2.*Znlk);
% end
% toc

%S2=zeros(size(S));
%tic;
denom=bsxfun(@times,(ht.^2).',Z);
nom1=bsxfun(@times,ht.',Z);
nom1=conj(nom1);
parfor k=1:K
    [lower_index,upper_index] = bsearch(Zis,k,k);
    if isempty(lower_index) || isempty(upper_index)
        continue;
    end
    ind=ZiI(lower_index:upper_index);
        nl=ceil(ind/KNN);
    Unlk=shift_phases(:,shifts(ind));  
    %Znlk=Z(ind);
    
    % XXX Why is ht of length N*L and not of length N?
    % XXX Why do you have imagninary valyes in the denominator?
    S(:,k)=sum(bsxfun(@times,(nom1(ind)).',Unlk.*X(:,nl)),2)/sum(denom(ind));
end
%toc
%aaa=1;

function [S,h,PCs,init_pick]=init(X,N,L,K)
% project all Fourier lines to the first nPCS principal componenets, 
% and then pick K signals (randomly) 
random_line=randi(L,[1,K]);
random_image=randi(N,[1,K-N]);
random_image=[(1:N) random_image];
init_pick= L*(random_image-1)+random_line;
init_pick=sort(init_pick(:));
nPCs=10;
%Xzm=X;
%m=mean(X,2);
%NL=size(X,2);
h=zeros(N*L,1);
%for k=1:NL
%    Xzm(:,k)=X(:,k)-m;
%end
%Xzm=cryo_raynormalize(Xzm);
[U,~]=svd(X);
PCs=U(:,1:nPCs);
S=PCs*PCs'*X;
S=cryo_raynormalize(S);
parfor n=1:N*L
    h(n)=real(S(:,n)'*X(:,n)); 
end
S=S(:,init_pick);
% %
% % Syntax:  [S,Z,Zi,PCs] = cld_ml(X,L,K,max_itr,noise_psd,verbose)
% %
% % Inputs:
% %    X - matrix of noisy polar fourier lines (vectors).
% %    L - angular resolution.
% %    K - number of clusters.
% %    KNN - number of clusters for the Estep (KNN<=K).
% %    max_itr - maximum number of iterations (maybe to be substitute with 
% %              a stoping condition...).
% %    noise_psd - noise power spectrum matrix.
% %    max_shift - maximum shift (in pixels).
% %    shift_step - shift step (in pixels).
% %
% % Outputs:
% %    S - clusters centers.
% %    Z - likelihoods.
% %    Zi - cluster indicators.
% %    PCs - principal componenets of X.
% %    mse - MSE in S (for debug)
% %
% % Example: 
% %
% % Other m-files required: none
% % Subfunctions: none
% % MAT-files required: none
% %
% % See also: 
% 
% % Author: Mor Cohen 
% 
% 
% narginchk(5, 9);
% mse_check=0;
% switch nargin
%     case 5
%         noise_psd=1;
%         verbose=1;
%         max_shift=0;
%         shift_step=1;
%     case 6
%         verbose=1;
%         max_shift=0;
%         shift_step=1;
%     case 7
%         max_shift=0;
%         shift_step=1;
%     case 8
%         error('shift parameter is missing');
%     case 10
%         mse_check=1;        
%         if verbose
%             fprintf('debug mode\n');
%         end
% end
% 
% if ~isinteger(K)
%     K=ceil(K);
% end
% 
% N=size(X,2)/L;
% if K<N
%     error('K input must be greater or equal N.')
% end
% NL=N*L;
% T=floor(max_shift/shift_step);
% shift_idx=(-T:T)*shift_step;
% N_shifts=length(shift_idx);
% C = inv(noise_psd);
% %% intialization (using PCA)
% Z=0;
% Zi=0;
% [S,h,shifts,PCs,init_pick]=init(X,N,L,K);
% mse=zeros(1,max_itr+2);
% if mse_check
%     Sclean=Xclean(:,init_pick);
%     Sclean=cryo_raynormalize(Sclean);
%     Y=cryo_raynormalize(X(:,init_pick));
%     d=abs(Y-Sclean).^2;
%     mse(1)=mean(d(:));
%     d=abs(S-Sclean).^2;
%     mse(2)=mean(d(:));
%     if verbose
%         msg=sprintf(' mse=%f\n',mse(2));
%         fprintf('%s',msg);
%     end
% end
% %% main loop
% for itr=1:max_itr
%     t1=clock;
%     if verbose
%         msg=sprintf('itr=%2d/%2d ',itr,max_itr);
%         fprintf('%s',msg);
%     end
%     
%     if verbose
%         fprintf('E step');
%     end
% 
%     Z=zeros(length(shift_idx),KNN,NL);
%     Zi=zeros(length(shift_idx),KNN,NL);
%     for n=1:N_shifts
%         shifts=repmat(shift_idx(n),NL,1);
%         [Z(n,:,:),Zi(n,:,:)]=Estep(sqrt(C)*X,L,shifts,h,sqrt(C)*S,KNN);
%     end
%     % taking just the top KNN elements in Z (out of n_shifts*KNN)
%     Z=reshape(Z, N_shifts*KNN,NL);
%     Zi=reshape(Zi, N_shifts*KNN,NL);
%     [~,ii]=sort(Z,1,'descend');
%     Zi_sorted=zeros(KNN,NL);
%     Z_sorted=zeros(KNN,NL);
%     for j = 1:size(Z,2)
%         Zi_sorted(:,j)=Zi(ii(1:KNN,j),j); 
%         Z_sorted(:,j)=Z(ii(1:KNN,j),j); 
%     end
%     Z=Z_sorted;
%     Zi=Zi_sorted;
%     clear Z_sorted;
%     clear Zi_sorted;
%     
%     if verbose
%         fprintf('/M step ');
%     end 
%     
%     shifts=shift_idx(ceil(ii(1:KNN,:)/KNN)); % Mstep_shifts
%     S=Mstep_S(X,Z,Zi,S,shifts,h,K,KNN);
%     h=Mstep_h(X,Z,Zi,S,shifts,h,C,N,L);
%     t2=clock;
%     t=etime(t2,t1);
%     if verbose
%         msg=sprintf(' t=%7.5f\n',t);
%         fprintf('%s',msg);
%     end
%     if mse_check
%         S=cryo_raynormalize(S);
%         d=abs(S-Sclean).^2;
%         mse(2+itr)=mean(d(:));
%         if verbose
%             msg=sprintf(' mse=%f\n',mse(2+itr));
%             fprintf('%s',msg);
%         end
%     end
% end
% if max_itr>0
%     Z=Z(1,:);
%     Zi=Zi(1,:);
% end
% 
% function [Z,Zi]=Estep(X,L,shiftst,ht,St,KNN)
% 
% NL=size(X,2);
% p=size(X,1);
% rho=(0:(p-1));
% Zi_cell=zeros(L,NL/L,KNN); % sliced variable
% Z_cell=zeros(L,NL/L,KNN); % sliced variable
% parfor n=1:NL/L
%     nl=(n-1)*L+(1:L);
%     shift_phases=exp(2*pi*sqrt(-1).*bsxfun(@times,rho',shiftst(nl)')./(2*p+1)); 
%     [Zi_cell(:,n,:),Z_cell(:,n,:)]=knn(bsxfun(@times,shift_phases.*X(:,nl),1./ht(nl)'),St,KNN);
% end
% Zi=reshape(Zi_cell,NL,KNN)';
% Z=reshape(Z_cell,NL,KNN);
% Z=bsxfun(@times,Z,ht.^2); % multiply by h^2 
% Z=exp(-Z');
% 
% function [h]=Mstep_h(X,Z,Zi,St,shiftst,h,C,N,L)
% p=size(X,1);
% rho=(0:(p-1));
% parfor n=1:N*L
%     Sk=St(:,Zi(:,n));
%     Znlk=Z(:,n);
%     shftnlk=shiftst(:,n);
%     shift_phases=exp(2*pi*sqrt(-1).*bsxfun(@times,rho',shftnlk')./(2*p+1)); 
%     Xnl=bsxfun(@times,shift_phases,X(:,n));
%     h(n)=sum(Znlk.*diag(real((C*Xnl)'*Sk)))/sum(Znlk.*diag(Sk'*(C*Sk)));
% end    
% 
% function S=Mstep_S(X,Z,Zi,S,shiftst,ht,K,KNN)
% p=size(X,1);
% rho=(0:(p-1));
% parfor k=1:K
%     ind=find(Zi==k);
%     if isempty(ind)
%         continue;
%     end
%     nl=ceil(ind/KNN);
%     Znlk=Z(ind);
%     sftnlk=shiftst(ind);
%     shift_phases=exp(2*pi*sqrt(-1).*bsxfun(@times,rho',sftnlk')./(2*p+1)); 
%     S(:,k)=sum(bsxfun(@times,(ht(nl).*Znlk)',shift_phases.*X(:,nl)),2)/sum(ht(nl).^2.*Znlk);
% end
% S=cryo_raynormalize(S);
% 
% function [S,h,shifts,PCs,init_pick]=init(X,N,L,K)
% % project all Fourier lines to the first nPCS principal componenets, 
% % and then pick K signals (randomly) 
% random_line=randi(L,[1,K]);
% random_image=randi(N,[1,K-N]);
% random_image=[(1:N) random_image];
% init_pick= L*(random_image-1)+random_line;
% init_pick=sort(init_pick(:));
% nPCs=10;
% Xzm=X;
% m=mean(X,2);
% NL=size(X,2);
% h=zeros(N*L,1);
% for k=1:NL
%     Xzm(:,k)=X(:,k)-m;
% end
% Xzm=cryo_raynormalize(Xzm);
% [U,~]=svd(Xzm);
% PCs=U(:,1:nPCs);
% S=PCs*PCs'*X;
% S=cryo_raynormalize(S);
% parfor n=1:N*L
%     h(n)=real(S(:,n)'*X(:,n)); 
% end
% S=S(:,init_pick);
% shifts=0;
% 
% 
