function [nThetaOpt,W_eps_opt ] = findOptEpsNTheta(x,ang_freqs,W_eps_init,nThetaThr)
%% Find optimal nTheta
N = size(x,2);
W = zeros(N,N,max(ang_freqs)+1);

% Correct distances to accomodate for negative frequencies
x_c = x;
x_c(ang_freqs~=0) = sqrt(2)*x_c(ang_freqs~=0);

%% Evaluate affinity matrices 
score = [];
i=1;
nTheta=32;
for m=0:max(ang_freqs)
    W(:,:,m+1) = x_c(ang_freqs==m,:)' * x_c(ang_freqs==m,:); 
end
Nnorm = bsxfun(@plus,sum(abs(x_c.').^2,2),sum(abs(x_c.').^2,2).');

WW = zeros(N);
while(1)
    WW_prev = WW;
    WW = exp(-( bsxfun(@minus,Nnorm,2*real(fft(W,nTheta,3))) )/W_eps_init);
    WW(repmat(logical(eye(size(WW,1))),1,1,size(WW,3))) = 0;    
    
    score(i) = max(max(abs(mean(WW,3)-mean(WW_prev,3))./abs(mean(WW,3))));
    if (score(i)<nThetaThr)
        break;
    end
    i=i+1;
    nTheta = nTheta+16;
end

nThetaOpt = nTheta + max(ang_freqs);
W_eps_opt = W_eps_init;

%% Find optimal eps
% lsr = 16;
% gridP = 30;
% W_range = logspace(log10(W_eps_init/lsr),log10(lsr*W_eps_init),gridP);
% 
% N = size(x,2);
% W = zeros(N,N,max(ang_freqs)+1);
% 
% % Correct distances to accomodate for negative frequencies
% x_c = x;
% x_c(ang_freqs~=0) = sqrt(2)*x_c(ang_freqs~=0);
% 
% %% Evaluate affinity matrices 
% score = zeros(1,gridP);
% for i=1:gridP
%     clc; i/gridP
%     for m=0:max(ang_freqs)
%         W(:,:,m+1) = x_c(ang_freqs==m,:)' * x_c(ang_freqs==m,:); 
%     end
% 
%     Nnorm = bsxfun(@plus,sum(abs(x_c.').^2,2),sum(abs(x_c.').^2,2).');
%     W = exp(-( bsxfun(@minus,Nnorm,2*real(fft(W,nTheta,3))) )/W_range(i));
%     W(repmat(logical(eye(size(W,1))),1,1,size(W,3))) = 0;    
%     
%     score(i) = log(sum(sum(sum(W))));
% end
% 
% secDiff = diff(diff(score));
% figure; plot(secDiff);

end

