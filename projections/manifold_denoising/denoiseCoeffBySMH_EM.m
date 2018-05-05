function [ coeff_denoise_curr,L,lScore] = denoiseCoeffBySMH_EM( y,x,ang_freqs,k,vCell,nv,maxIter,nTheta )
if (numel(k)==1)
    k = repmat(k,1,max(ang_freqs)+1);
end

fftLen = max(ang_freqs)+1;
nTheta = max(nTheta,fftLen);
coeff_denoise_curr = x;
Z0 = rand(1,size(x,2));

L=zeros(1,maxIter);
chunkSize = 128;
nChunks = ceil(size(y,2)/chunkSize);
mIdx = unique(ang_freqs).';

for itrCnt=1:maxIter    
%     itrCnt
    disp(['Iteration: ',num2str(itrCnt)]);
    Lcurr = 0;
    lScore = zeros(1,size(x,2));
    lScore2 = zeros(1,size(x,2));
    A = cell(1,max(ang_freqs)+1);
    B = cell(1,max(ang_freqs)+1);
    VZ = cell(nChunks,max(ang_freqs)+1);
%     VZ = cell(1,max(Freqs)+1);
    Z0 = zeros(1,size(x,2));
%     VZ = zeros(size(pCoeff(Freqs==N,:),1),size(coeff_denoise_curr,2));
    for m=0:max(ang_freqs)
            freqIdx = m+1;
            A{freqIdx} = zeros(k(freqIdx));
            B{freqIdx} = zeros(size(y(ang_freqs==m,:),1),k(freqIdx));
            for i=1:nChunks
                VZ{i,freqIdx} = zeros(size(y(ang_freqs==m,:),1),size(coeff_denoise_curr,2));
            end
%             VZ{freqIdx} = zeros(size(pCoeff(Freqs==N,:),1),size(coeff_denoise_curr,2));
    end
    parfor i = 1:nChunks
        currIdx = ((i-1)*chunkSize+1):min(i*chunkSize,size(y,2));    
%         currIdx(end)
        Zcurr = zeros(numel(currIdx),size(coeff_denoise_curr,2),max(ang_freqs)+1);
        for m=mIdx
%             clc; j/max(Freqs)
            Zcurr(:,:,m+1) = y(ang_freqs==m,currIdx)' * coeff_denoise_curr(ang_freqs==m,:); 
        end
        Nnorm = (bsxfun(@plus,sum(abs(y(:,currIdx).').^2,2),sum(abs(coeff_denoise_curr.').^2,2).'));
        Zcurr = bsxfun(@minus,Nnorm,2*real(nTheta*ifft(Zcurr,nTheta,3)));
        minDcurr = min(min(Zcurr,[],3),[],2);
        %Zcurr_score = exp(-Zcurr/nv/2 + size(y,1)/2);
        Zcurr = exp(-(Zcurr-minDcurr)/nv/2);
%         Zcurr(1:(chunkSize+1):chunkSize^2) = 0;
        rowSum = sum(sum(Zcurr,3),2);
        for j=1:nTheta
            Zcurr(:,:,j) = bsxfun(@times,Zcurr(:,:,j),1./rowSum); 
        end
        lScore = lScore + sum(sum(Zcurr,3),1);        
%         lScore = max([max(sum(Zcurr,3),[],1); lScore],[],1);
        Zcurr = fft(Zcurr,nTheta,3);        
        Zcurr = Zcurr(:,:,1:fftLen);
        
        Lcurr = Lcurr + sum(log(rowSum)-minDcurr/nv/2);
        
        Z0 = Z0 + sum(Zcurr(:,:,1));
        tmp = cell(1,max(ang_freqs)+1);
        for m=mIdx
            freqIdx = m+1;            
            tmp{freqIdx} = y(ang_freqs==m,currIdx) * Zcurr(:,:,freqIdx);
        end
        VZ(i,:) = tmp;
    end
    L(itrCnt) = Lcurr;
    
    VZtot = cell(1,max(ang_freqs)+1);
    for m=mIdx
        freqIdx = m+1;
        VZtot{freqIdx} = zeros(size(y(ang_freqs==m,:),1),size(coeff_denoise_curr,2));
    end
    for m=mIdx
        freqIdx = m+1;
        for i=1:nChunks
            VZtot{freqIdx} = VZtot{freqIdx} + VZ{i,freqIdx};
        end
    end        
    
    for m=mIdx
%         m
        freqIdx = m+1;
        A{freqIdx} = vCell{freqIdx}(1:size(x,2),1:k(freqIdx)).' * bsxfun(@times,Z0.',conj(vCell{freqIdx}(1:size(x,2),1:k(freqIdx))));
        B{freqIdx} = VZtot{freqIdx} * conj(vCell{freqIdx}(1:size(x,2),1:k(freqIdx)));                          
        coeffCell{m+1} = B{freqIdx}/A{freqIdx};  
        coeff_denoise_curr(ang_freqs==m,:) = coeffCell{freqIdx}(:,1:k(freqIdx)) * vCell{freqIdx}(1:size(x,2),1:k(freqIdx)).';        
    end       
        
end

end

