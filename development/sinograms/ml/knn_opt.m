function [N,D] = knn_opt(X, Y, k)
% [ index,distance] = knn(query, data, k);
% top(x,k) is a partial sort function which returns top k entries of x %
%[D, N] = top(bsxfun(@plus,(-2)*real(Y'*X),dot(Y,Y,1)'), k);

YtY=sum(abs(Y).^2);
YtY=YtY.';
D1=bsxfun(@plus,(-2)*real(Y'*X),YtY)';
% get the top k elements from each column
[~,idx]=sort(D1,2);
N=idx(:,1:k);
N=N';
D=zeros(size(D1,1),k);
for j=1:size(D1,1)
    D(j,:)=D1(j,idx(j,1:k));
end
D=D';
XtX=sum(abs(X).^2);
D = bsxfun(@plus,D,XtX);
