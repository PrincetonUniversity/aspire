function [N,D] = knn_opt_gpu_lowmem(X, Y, knn,MAXMEM)
% [ index,distance] = knn(query, data, k);
% top(x,k) is a partial sort function which returns top k entries of x %
%[D, N] = top(bsxfun(@plus,(-2)*real(Y'*X),dot(Y,Y,1)'), k);


% %%%% Reference code
% Y=single(Y);
% X=single(X);
% YtY_ref=sum(abs(Y).^2);
% YtY_ref=YtY_ref.';
% D1_ref=bsxfun(@plus,(-2)*real(Y'*X),YtY_ref)';
% % get the top k elements from each column
% [~,idx]=sort(D1_ref,2);
% N_ref=idx(:,1:knn);
% N_ref=N_ref';
% D_ref_1=zeros(size(D1_ref,1),knn);
% for j=1:size(D1_ref,1)
%     D_ref_1(j,:)=D1_ref(j,idx(j,1:knn));
% end
% D_ref_1=D_ref_1';
% XtX_ref=sum(abs(X).^2);
% D_ref = bsxfun(@plus,D_ref_1,XtX_ref);
% %%%% Rerference code

if nargin<4
    MAXMEM=1.0e9; % 1GB on the GPU
end

% We need to solve for blocksize that satisfies 
%   (numel(X)+size(Y,1)*blocksize+blocksize*size(X,2))*4<MAXMEM

numbytes=4; % Single precision
blocksize=floor((MAXMEM/numbytes-numel(X))/(size(Y,1)+size(X,2)));

if blocksize<100
    error('Insufficient memory. Allocate more memory');
end

D1=zeros(size(X,2),size(Y,2),'single'); % Allocate memory
gX=gpuArray(X);
Niter=ceil(size(Y,2)/blocksize);

for k=1:Niter
    idx=(k-1)*blocksize+1:min(k*blocksize,size(Y,2));
    gY=gpuArray(Y(:,idx));

    gYtY=sum(abs(gY).^2);
    gYtY=gYtY.';
    gD1=bsxfun(@plus,(-2)*real(gY'*gX),gYtY)';
    D1(:,idx)=gather(gD1);
end

% get the top k elements from each column

% The following three lines were replaced by the line following it.
% [~,idx]=sort(D1,2);
%  N=idx(:,1:knn);
%  N=N';

[~,N]=mink(double(D1.'),knn);


 
 D=zeros(size(D1,1),knn);
 for j=1:size(D1,1)
     D(j,:)=D1(j,N(1:knn,j)); 
 end
 D=D';
 XtX=sum(abs(X).^2);
 D = bsxfun(@plus,D,XtX);
