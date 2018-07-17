function [vCell,dCell,denMat] = evalSMH(x,ang_freqs,W_eps,nTheta,maxEigIdx,t_Thr,normalizeDensity)
N = size(x,2);
W = zeros(N,N,max(ang_freqs)+1);
fftLen = max(ang_freqs)+1;

% Correct distances to accomodate for negative frequencies
x_c = x;
x_c(ang_freqs~=0) = sqrt(2)*x_c(ang_freqs~=0);

%% Evaluate affinity matrics 
for i=0:max(ang_freqs)
%     i
    W(:,:,i+1) = x_c(ang_freqs==i,:)' * x_c(ang_freqs==i,:); 
end

Nnorm = bsxfun(@plus,sum(abs(x_c.').^2,2),sum(abs(x_c.').^2,2).');
W = exp(-( bsxfun(@minus,Nnorm,2*real(fft(W,nTheta,3))) )/W_eps);
W(repmat(logical(eye(size(W,1))),1,1,size(W,3))) = 0;

%% Compute the Fourier transform of W
W = fft(W,nTheta,3);
W = W(:,:,1:fftLen);
W0 = W(:,:,1);
D = diag(1./sum(W0,2));
denMat = D;

%% Normalize density on manifold
if normalizeDensity    
    for i=0:max(ang_freqs)
        W(:,:,i+1) = W(:,:,i+1)*denMat;
    end
    W0 = W(:,:,1);
    D = diag(1./sum(W0,2));
end

%% Construct graph laplacian and find eigenvectors
vCell = cell(1,max(ang_freqs)+1);
dCell = cell(1,max(ang_freqs)+1);
parfor i=0:max(ang_freqs)
%     i    
    if ~isempty(maxEigIdx)    
        P = D*W(:,:,i+1);
        P(abs(P)<t_Thr*abs(max(P,[],2))*ones(1,N)) = 0;
        [v,d] = eigs(sparse(P),maxEigIdx,'lr');
    else        
        P = sqrt(D)*W(:,:,i+1)*sqrt(D);
        P = (P+P')/2;   % Enforce Hermitian matrix numerically
        [v,d] = eig(P);
        v = sqrt(D)*v;
    end
    [d,idx] = sort(diag(real(1-d)),'ascend');
    v=v(:,idx);
    vCell{i+1} = v;
    dCell{i+1} = d;
end


end

