function [NN,D] = knn_opt_gpu_lowmem_2(X, Y, knn,MAXMEM)
% [ index,distance] = knn(query, data, k);
% top(x,k) is a partial sort function which returns top k entries of x %
%[D, N] = top(bsxfun(@plus,(-2)*real(Y'*X),dot(Y,Y,1)'), k);

if nargin<4
    MAXMEM=1.0e9; % 1GB on the GPU
end

% We need to solve for blocksize that satisfies 
%   (numel(X)+size(Y,1)*blocksize+blocksize*size(X,2))*4<MAXMEM

numbytes=4; % Single precision
blocksize=floor((MAXMEM/numbytes-numel(X))/(size(Y,1)+size(X,2))/2);

if blocksize<100
    error('Insufficient memory. Allocate more memory');
end

% DEBUG
% D1_debug=zeros(size(X,2),size(Y,2),'single'); % Allocate memory
% DEBUG
gX=gpuArray(single(X));
Niter=ceil(size(Y,2)/blocksize);
gcandNN=gpuArray(single(zeros(knn*Niter,size(X,2)))); % Candidate nearest neighbors
gcandD=gpuArray(single(zeros(size(gcandNN))));        % Distance  to candidate NN
for k=1:Niter
    idx=(k-1)*blocksize+1:min(k*blocksize,size(Y,2));
    gY=gpuArray(single(Y(:,idx)));

    gYtY=sum(abs(gY).^2);
    gYtY=gYtY.';
    gD1=bsxfun(@plus,(-2)*real(gY'*gX),gYtY)';
    
    % DEBUG
    % D1_debug(:,idx)=gather(gD1);
    % DEBUG
    %tic;
    gD1=gD1.';
    maxval=max(gD1(:))+1; % Use this value to reset each candidate nearest 
                          % neighbor we find, so we won't find it again.

    % Find knn nearest neighbors candidates using the current block point
    % Y. Later we will find the true nearest neighbors among these
    % candidates. 
    for j=1:knn
        [gCandVal,gCandIdx]=min(gD1); % Find the next nn cadidate in each column.
        %gidx=sub2ind(size(gD1),gCandIdx,1:size(X,2)); % convert the index of the 
        gidx=gCandIdx+(0:size(X,2)-1).*size(gD1,1); % Fast implementation of the above row.
            % current candidate to a linear index.
        gD1(gidx)=maxval;   % Reset the distance of all candidates 
                            % detected, so we don't detect again in the
                            % next iteration.

        gcandNN((k-1)*knn+j,:)=gCandIdx+(k-1)*blocksize; 
        gcandD((k-1)*knn+j,:)=gCandVal; 
    end
    %toc    
end

% get the top knn elements from each column
[gSortedCandD,gSortedidx]=sort(gcandD);
Sortedidx=gather(gSortedidx);
candNN=gather(gcandNN);
NN=zeros(knn,size(X,2));
for k=1:size(X,2)
    NN(:,k)=candNN(Sortedidx(1:knn,k),k);    
end

% Calculate distance to the nearest neighbors
D=gather(gSortedCandD);
D=D(1:knn,:);
XtX=sum(abs(X).^2);
D = bsxfun(@plus,D,XtX);

% DEBUG
%
% [~,N_debug]=mink(D1_debug.',knn);
% D_debug=zeros(size(D1_debug,1),knn);
% for j=1:size(D1_debug,1)
%     D_debug(j,:)=D1_debug(j,N_debug(:,j));
% end
% D_debug=D_debug.';
% XtX=sum(abs(X).^2);
% D_debug = bsxfun(@plus,D_debug,XtX);
% 
% disp(norm(D(:)-D_debug(:))./norm(D_debug(:))) % This should be zero
% disp(norm(NN(:)-N_debug(:))/norm(N_debug(:))) % This may not be zero are
% % there are elements with numerical multiplicity in D and so NN and
% % NN_debug may use different indices to point a given distance in D. In
% % any case, the value of this norm should be tiny (10^-7).
% DEBUG