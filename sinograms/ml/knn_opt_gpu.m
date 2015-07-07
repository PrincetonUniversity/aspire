function [N,D] = knn_opt_gpu(X, Y, k)
% [ index,distance] = knn(query, data, k);
% top(x,k) is a partial sort function which returns top k entries of x %
%[D, N] = top(bsxfun(@plus,(-2)*real(Y'*X),dot(Y,Y,1)'), k);

gX=gpuArray(single(X));
gY=gpuArray(single(Y));

gYtY=sum(abs(gY).^2);
gYtY=gYtY.';
gD1=bsxfun(@plus,(-2)*real(gY'*gX),gYtY)';
% get the top k elements from each column
D1=gather(gD1);
% The following three lines were replaced by the line following it.
% [~,idx]=sort(D1,2);
%  N=idx(:,1:k);
%  N=N';

 [~,N]=mink(D1.',k);
 D=zeros(size(D1,1),k);
 for j=1:size(D1,1)
     D(j,:)=D1(j,N(:,j)); 
 end
 D=D';
 XtX=sum(abs(X).^2);
 D = bsxfun(@plus,D,XtX);
