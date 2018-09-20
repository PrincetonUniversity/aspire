function [vCell,dCell] = evalSMH_sparse(x,ang_freqs,W_eps,nTheta,maxEigIdx,nn,dispFlag)
N = size(x,2);
maxEigIdx = min(maxEigIdx,N);
% W = zeros(N,N,max(ang_freqs)+1);
fftLen = max(ang_freqs)+1;
mIdx = (unique(ang_freqs)+1).';

% Correct distances to accomodate for negative frequencies
% x_c = (x);
x_c = gpuArray(single(x));
x_c(ang_freqs~=0) = sqrt(2)*x_c(ang_freqs~=0);

%% Compute steerable Graph-Laplacian quantities
chunkSize = 64;
nChunks = ceil(N/chunkSize);

W_nn = gpuArray(single(zeros(N,nn,fftLen)));
W_nn_idx = gpuArray(zeros(N,nn));

normVec = sum(abs(x_c.').^2,2);
for j = 1:nChunks
%     j/nChunks
    if dispFlag
        %clc; 
	%disp('Computing steerable graph Laplacian pair-wise affinities...');
        log_message(['Computing chunk ',num2str(j),' out of ',num2str(nChunks)]);
    end
    currIdx = ((j-1)*chunkSize+1):min(j*chunkSize,N); 
%     W = zeros(numel(currIdx),N,max(ang_freqs)+1);
    W = gpuArray(single(zeros(numel(currIdx),N,max(ang_freqs)+1)));
    W_nn_curr = gpuArray(single(zeros(numel(currIdx),nn,nTheta)));
    %% Evaluate affinity matrices 
    for i=(mIdx-1)
        W(:,:,i+1) = x_c(ang_freqs==i,currIdx)' * x_c(ang_freqs==i,:); 
    end
    Nnorm = bsxfun(@plus,normVec(currIdx),normVec.');
    W = exp(-( bsxfun(@minus,Nnorm,2*real(fft(W,nTheta,3))) )/W_eps);
    diagInd = logical(eye(N));
    diagInd = repmat(diagInd(currIdx,:),1,1,size(W,3));
    W(diagInd) = 0;
    %% Construct nearest-neighbour weight matrix
    W0_curr = sum(W,3);
    [~,sortIdx] = sort(W0_curr.','descend'); sortIdx=sortIdx.';
    for i=1:numel(currIdx)
        W_nn_curr(i,:,:) = W(i,sortIdx(i,1:nn),:);
%         W_nn(currIdx(i),:,:) = W(i,sortIdx(i,1:nn),:);
    end
    W_nn_idx(currIdx,:) = sortIdx(:,1:nn);
    %% Compute the Fourier transform of W
%     W = fft(W,nTheta,3);
%     W = W(:,:,1:fftLen);
%     W0_curr = real(W(:,:,1));
    W_nn_curr = fft(W_nn_curr,nTheta,3);
    W_nn(currIdx,:,:) = W_nn_curr(:,:,1:fftLen);
    clear W;
end

%% Construct graph laplacian and find eigenvectors
vCell = cell(1,max(ang_freqs)+1);
dCell = cell(1,max(ang_freqs)+1);
% nnp = (nnp+nnp.')>0;
% nnp = nnp.*(nnp.');
% W0 = real(outOfMemFile{1});
% D = 1./sum(W0.*nnp,2);
W0_nn = real(W_nn(:,:,1));
W_nn_curr = double(W0_nn);
W0 = gather(sparse(repmat((1:N).',nn,1),W_nn_idx(:),W_nn_curr(:),N,N));
% D = 1./sum(W0_nn,2);
pattern = (gather(W0)>0);
pattern_sym = ((pattern-pattern.')~=0);
W0_sym = (W0+W0')/2;
W0_sym(pattern_sym) = W0_sym(pattern_sym)*2;
D = gpuArray(1./sum(W0_sym,2));
for i=mIdx
    if dispFlag
        %clc; 
	%disp('Evaluating the steerable manifold harmonics...');
        log_message(['Evaluating angular index ',num2str(i),' out of ',num2str(max(mIdx))]);
    end
%     i
%     P = outOfMemFile{i+1};
%     P = P.*nnp;     % Enforce sparseness
    W_nn_curr = double(W_nn(:,:,i));
    P = sparse(repmat((1:N).',nn,1),W_nn_idx(:),W_nn_curr(:),N,N);
%     P = bsxfun(@times,P,sqrt(D));
%     P = bsxfun(@times,P.',sqrt(D)).';
    P = gather(P);
    P = (P+P')/2;
    P(pattern_sym) = P(pattern_sym)*2;
    P = gpuArray(P);
    D_s = spdiags(sqrt(D),0,N,N);
    P = D_s*P*D_s;
%     P = speye(N)-D_s*P*D_s;
%     P = sqrt(D)*W(:,:,i+1)*sqrt(D);
%     P = gpuArray((P+P')/2);   % Enforce Hermitian matrix numerically  
%     P = gpuArray((P+P')/2);   % Enforce Hermitian matrix numerically 
    P=single(full(P));
    [v,d] = eig((P+P')/2);
%     [v,s,~] = svds(2*speye(N)-P,maxEigIdx);
%     d = 2-s;
%     v = sqrt(D)*v;
    v = bsxfun(@times,gather(v),sqrt(D));
    [d,sortIdx] = sort(gather(real(1-diag(d))),'ascend');
    v=v(:,sortIdx);
    vCell{i} = gather(v(:,1:maxEigIdx));
    dCell{i} = gather(d(1:maxEigIdx));
%     vCell{i+1} = v;
%     dCell{i+1} = diag(d);
    clear P;
end


end

