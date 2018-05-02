function [vCell,dCell] = evalSMHOutOfMem_v2(x,ang_freqs,W_eps,nTheta,maxEigIdx,t_Thr)
N = size(x,2);
% W = zeros(N,N,max(ang_freqs)+1);
fftLen = max(ang_freqs)+1;

% Correct distances to accomodate for negative frequencies
x_c = (x);
x_c(ang_freqs~=0) = sqrt(2)*x_c(ang_freqs~=0);

%% Compute steerable Graph-Laplacian quantities
chunkSize = 128;
nChunks = ceil(N/chunkSize);

W = single(zeros(N)+1i*eps*ones(N));
outOfMemFile = cell(fftLen,1);
for i=1:1:fftLen
    currFile = ['./rotInvImManDeNoise/outOfMemFiles/sGLquan_',num2str(i)];
    save(currFile,'W','-v7.3');
    outOfMemFile{i} = matfile(currFile,'Writable',true);
end

W = single(zeros(chunkSize,N,fftLen)+1i*eps*ones(chunkSize,N,fftLen));
D = single(zeros(N,1));
normVec = sum(abs(x_c.').^2,2);
for j = 1:nChunks
    currIdx = ((j-1)*chunkSize+1):min(j*chunkSize,N); 
    W = zeros(numel(currIdx),N,max(ang_freqs)+1);
    %% Evaluate affinity matrics 
    for i=0:max(ang_freqs)
        W(:,:,i+1) = x_c(ang_freqs==i,currIdx)' * x_c(ang_freqs==i,:); 
    end
    Nnorm = bsxfun(@plus,normVec(currIdx),normVec.');
    W = exp(-( bsxfun(@minus,Nnorm,2*real(fft(W,nTheta,3))) )/W_eps);
    diagInd = logical(eye(N));
    diagInd = repmat(diagInd(currIdx,:),1,1,size(W,3));
    W(diagInd) = 0;
    %% Compute the Fourier transform of W
    W = fft(W,nTheta,3);
    W = W(:,:,1:fftLen);
    W0 = real(W(:,:,1));
    D(currIdx) = 1./sum(W0,2);
    %% Append quantities to files
    for i=1:1:fftLen
        outOfMemFile{i}.W(currIdx,:) = single(W(:,:,i));
    end
end

%% Construct graph laplacian and find eigenvectors
vCell = cell(1,max(ang_freqs)+1);
dCell = cell(1,max(ang_freqs)+1);
for i=0:max(ang_freqs)
    i
    P = outOfMemFile{i+1}.W;
    P = bsxfun(@times,P,sqrt(D));
    P = bsxfun(@times,P.',sqrt(D)).';
%     P = sqrt(D)*W(:,:,i+1)*sqrt(D);
    P = gpuArray((P+P')/2);   % Enforce Hermitian matrix numerically    
    [v,d] = eig(P);
%     v = sqrt(D)*v;
    v = bsxfun(@times,gather(v),sqrt(D));
    [d,sortIdx] = sort(gather(real(1-diag(d))),'ascend');
    v=v(:,sortIdx);
    vCell{i+1} = v(:,1:maxEigIdx);
    dCell{i+1} = d(1:maxEigIdx);
end


end

