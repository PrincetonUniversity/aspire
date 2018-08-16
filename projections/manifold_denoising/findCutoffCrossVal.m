function lambdaOpt = findCutoffCrossVal(x,ang_freqs,vCell,dCell,dl,nTheta,nv,trainSize)
%% Find frequency cutoff by empirical likelihood cross-validation 
lambda_c = dl:dl:1-dl;
% lambda_c = 0.7:dl:1-dl;
empLkh = -inf(1,numel(lambda_c));

%clc; 
disp('Performing cross-validation for optimal cut-off frequency...')
for k=1:numel(lambda_c)
    empLkh(k) = getCVlogL_gpu( x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c(k));
    disp(['Log-likelihood for cut-off frequency ',num2str(lambda_c(k)),': ',num2str(empLkh(k))]);
%     if (k>1)
%         if empLkh(k)<empLkh(k-1)    % Break loop if there is a decrease in the empirical likelihood
%             break;
%         end
%     end
end
    
% figure; plot(empLkh); grid on;

[~,kOpt] = max(empLkh);
lambdaOpt = lambda_c(kOpt);

end

