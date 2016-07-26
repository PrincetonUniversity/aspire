function [Z,Zi]=Estep(X,L,ht,St,KNN)

assert(~any(isnan(ht)))    % ht must not contain NaNs.
assert(~any(isnan(St(:)))) % St must not contain NaNs.
NL=size(X,2);
N=ceil(NL/L);
Zi=zeros(KNN,L,N); % sliced variable
Z=zeros(KNN,L,N); % sliced variable
for n=1:N
    %n
    nl=(n-1)*L+(1:L);
    [Zi(:,:,n),Z(:,:,n)]=knn_opt_gpu_low_mem_2(bsxfun(@times,X(:,nl),1./ht(nl)'),St,KNN,0.5e9);
    %[Zi(:,:,n),Z(:,:,n)]=knn_opt_gpu_lowmem(bsxfun(@times,X(:,nl),1./ht(nl)'),St,KNN,0.5e9);
    
    % COMAPRE TO PREVIOUS VERSION
    %tic
    %[A,B]=knn_opt_gpu_lowmem(bsxfun(@times,X(:,nl),1./ht(nl)'),St,KNN,0.5e9);   
    %[A,B]=knn_opt_gpu(bsxfun(@times,X(:,nl),1./ht(nl)'),St,KNN);
    %toc
    %disp(norm(Z(:,:,n)-B)/norm(Z(:,:,n))) % Should be zero
    %disp(norm(Zi(:,:,n)-A)/norm(Zi(:,:,n))) % Should be about 10^-7 for single precision
    % COMPARE TO PREVIOUS VERSION
    %[Zi(:,:,n),Z(:,:,n)]=knn_opt_gpu(bsxfun(@times,X(:,nl),1./ht(nl)'),St,KNN);
end

Zi=reshape(Zi,KNN,NL);
Z=reshape(Z,KNN,NL);
Z=bsxfun(@minus,Z,min(Z));
Z=bsxfun(@times,Z,(ht').^2); % multiply back by h^2 
Z=exp(-Z);
assert(norm(Z(:))>0);
