function [ vCell_ext ] = smhExtNystrom(x,vCell,dCell,nImSMH,ang_freqs,W_eps,nTheta,normalizeDensity,denMat)
%% Correct distances to accomodate for negative frequencies
x_c = x; clear x;
x_c(ang_freqs~=0) = sqrt(2)*x_c(ang_freqs~=0);

%% Extend eigenvectors to all points by Nystrom's method
mIdx = (unique(ang_freqs)).';
fftLen = max(ang_freqs)+1;
chunkSize = 128;
vCell_ext = vCell;
originalIdx = 1:nImSMH;
extendIdx = setdiff(1:size(x_c,2),originalIdx);
for j=0:max(ang_freqs)
    vCell_ext{j+1} = zeros(size(x_c,2),size(vCell{j+1},2));
    vCell_ext{j+1}(originalIdx,:) = vCell{j+1};
    vCell_ext{j+1}(extendIdx,:) = zeros(numel(extendIdx),size(vCell{j+1},2));
end 
nChunks = ceil(numel(extendIdx)/chunkSize);
x_c_origin = x_c(:,originalIdx);
x_c_extend = x_c(:,extendIdx);
clear x_c;
for i = 1:nChunks
    currIdx = ((i-1)*chunkSize+1):min(i*chunkSize,numel(extendIdx));    
    currIdx(end)
    Zcurr = zeros(numel(currIdx),numel(originalIdx),nTheta);
    for m=0:max(ang_freqs)
        Zcurr(:,:,m+1) = x_c_extend(ang_freqs==m,currIdx)' * x_c_origin(ang_freqs==m,:); 
    end
    Nnorm = (bsxfun(@plus,sum(abs(x_c_extend(:,currIdx).').^2,2),sum(abs(x_c_origin.').^2,2).'));
    Zcurr = exp(-( bsxfun(@minus,Nnorm,2*real(fft(Zcurr,nTheta,3))) )/W_eps);
    Zcurr = Zcurr + eps; % Prevents zero rows due to outliers.
    
    Zcurr = fft(Zcurr,nTheta,3);
    Zcurr = Zcurr(:,:,1:fftLen);
    
    if normalizeDensity
        Zcurr = bsxfun(@times,Zcurr,diag(denMat).');
    end
        
    D = diag(1./sum(Zcurr(:,:,1),2));
    for j=mIdx
        vCell_ext{j+1}(extendIdx(currIdx),:) = D*(Zcurr(:,:,j+1)*vCell{j+1}*diag(1./(1-dCell{j+1})));
    end
end

end

