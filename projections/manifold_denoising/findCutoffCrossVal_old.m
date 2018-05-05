function lambdaOpt = findCutoffCrossVal_old(x,ang_freqs,vCell,dCell,dl,nTheta,nv,trainSize)
%% Find frequency cutoff by empirical likelihood cross-validation 
% kVec = 1:dk:size(evIdxMat,1);
lambda_c = dl:dl:1-dl;
% chunkSize = 128;
empLkh = -inf(1,numel(lambda_c));

% nImages = size(x,2);
% % trainIdx = 1:ceil(nImages/2);   % Use half the data for train and half for testing
% trainIdx = 1:trainSize;
% testIdx = trainIdx(end)+1:nImages;
% y = x(:,testIdx); 
% nChunks = ceil(numel(testIdx)/chunkSize);
for k=1:numel(lambda_c)
%     clc; k/numel(kVec)
%     kCurr = find(lambda>lambda_c(k),1,'first')-1;
%     if kCurr>0
% %     [coeff_cv,rank] = denoiseCoeffBySMH(x,ang_freqs,evIdxMat(kVec(k),:),vCell_ext,nv,trainIdx);
%         [coeff_cv,rank] = denoiseCoeffBySMH(x,ang_freqs,evIdxMat(kCurr,:),vCell_ext,nv,trainIdx);
%         maxEmpF = find(rank>0,1,'last')-1;
%         empLkh(k) = 0;
%         for i = 1:nChunks
%             currIdx = ((i-1)*chunkSize+1):min(i*chunkSize,numel(testIdx));    
%             Zcurr = zeros(numel(currIdx),numel(trainIdx),maxEmpF+1);
%             for m=0:maxEmpF
%                 Zcurr(:,:,m+1) = y(ang_freqs==m,currIdx)' * coeff_cv(ang_freqs==m,:); 
%             end
%             Nnorm = (bsxfun(@plus,sum(abs(y(:,currIdx).').^2,2),sum(abs(coeff_cv.').^2,2).'));
%             Zcurr = bsxfun(@minus,Nnorm,2*real(nTheta*ifft(Zcurr,nTheta,3)));
%             Zcurr = Zcurr - size(y,1)*nv;   % Normalize for better numerics
%             Zcurr = exp(-Zcurr/nv/2);
%             empLkh(k) = empLkh(k) + sum(log(sum(sum(Zcurr,3),2)));
%     %         empLkh(k) = empLkh(k) + sum(log(sum(sum(exp(-bsxfun(@minus,Nnorm,2*real(nTheta*ifft(Zcurr,nTheta,3)))/nv/2),3),2)));
%         end
    empLkh(k) = getCVlogL_gpu( x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c(k));
    disp(['Log-likelihood for cut-off frequency',num2str(lambda_c(k)),': ',num2str(empLkh(k))]);
    if (k>1)
        if empLkh(k)<empLkh(k-1)    % Break loop if there is a decrease in the empirical likelihood
            break;
        end
    end
end
    
% figure; plot(empLkh); grid on;

[~,kOpt] = max(empLkh);
lambdaOpt = lambda_c(kOpt);

end

