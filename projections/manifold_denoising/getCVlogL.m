function [ empLkh ] = getCVlogL( x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c)
%% Sort the steerable manifold harmonics according to their frequencies
dMat = [];
for i=0:max(ang_freqs)
    dMat = [dMat dCell{i+1}];
end
freqMat = repmat(0:max(ang_freqs),numel(dCell{1}),1);
[lambda,dSortIdx] = sort(dMat(:),'ascend');
freqMat = freqMat(:); freqMat = freqMat(dSortIdx);

evIdxMat = zeros(numel(lambda),max(ang_freqs)+1);
evIdxMat(1,freqMat(1)+1) = 1;
for i=2:numel(freqMat)
    evIdxMat(i,:) = evIdxMat(i-1,:);
    evIdxMat(i,freqMat(i)+1) = evIdxMat(i,freqMat(i)+1) + 1;
end

%% Compute log-likelihood
chunkSize = 128;

nImages = size(x,2);
trainIdx = 1:trainSize;
testIdx = trainIdx(end)+1:nImages;
y = x(:,testIdx); 
nChunks = ceil(numel(testIdx)/chunkSize);
kCurr = find(lambda>lambda_c,1,'first')-1;
empLkh = -inf;
if kCurr>0
%     [coeff_cv,rank] = denoiseCoeffBySMH(x,ang_freqs,evIdxMat(kVec(k),:),vCell_ext,nv,trainIdx);
    [coeff_cv,rank] = denoiseCoeffBySMH(x,ang_freqs,evIdxMat(kCurr,:),vCell,nv,trainIdx);
    maxEmpF = find(rank>0,1,'last')-1;
    empLkh = 0;
    for i = 1:nChunks
        currIdx = ((i-1)*chunkSize+1):min(i*chunkSize,numel(testIdx));    
        Zcurr = zeros(numel(currIdx),numel(trainIdx),maxEmpF+1);
        for m=0:maxEmpF
            Zcurr(:,:,m+1) = y(ang_freqs==m,currIdx)' * coeff_cv(ang_freqs==m,:); 
        end
        Nnorm = (bsxfun(@plus,sum(abs(y(:,currIdx).').^2,2),sum(abs(coeff_cv.').^2,2).'));
        Zcurr = bsxfun(@minus,Nnorm,2*real(nTheta*ifft(Zcurr,nTheta,3)));
        Zcurr = Zcurr - size(y,1)*nv;   % Normalize for better numerics
        Zcurr = exp(-Zcurr/nv/2);
        empLkh = empLkh + sum(log(sum(sum(Zcurr,3),2)));
    end
end    


end

