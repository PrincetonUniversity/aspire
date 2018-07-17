function [x_hat,rank] = denoiseCoeffBySMH(x,ang_freqs,k,vCell,nv,dataValidIdx)
if (numel(k)==1)
    k = repmat(k,1,max(ang_freqs)+1);
end
coeffCell = cell(1,max(ang_freqs));
x_hat = zeros(numel(ang_freqs),numel(dataValidIdx));
% x_hat = zeros(numel(ang_freqs),size(x,2));
rank = zeros(1,max(ang_freqs)+1);
mIdx = unique(ang_freqs).';
for m=mIdx
    freqIdx = m+1;
    [q,~] = qr(vCell{freqIdx}(dataValidIdx,1:k(freqIdx)),0);
    coeffCell{freqIdx} = x(ang_freqs==m,dataValidIdx)*conj(q);
%     [coeffCell{freqIdx},rank(freqIdx)] = matrixDenoise(coeffCell{freqIdx},nv);
    rank(freqIdx) = min(size(q));
    x_hat(ang_freqs==m,:) = coeffCell{freqIdx}(:,1:k(freqIdx)) * q.';
%     x_hat(ang_freqs==m,:) = coeffCell{freqIdx}(:,1:k(freqIdx)) * q.' + matrixDenoise(x(ang_freqs==m,dataValidIdx)-(x(ang_freqs==m,dataValidIdx)*conj(q))*q.',nv);
end
end

