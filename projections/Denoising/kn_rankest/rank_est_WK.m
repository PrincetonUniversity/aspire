function K = rank_est_WK(lambda,p,n)
% function K = rank_est_WK(lambda,p,n)
% MDL estimator for number of signals proposed by
% Wax & Kailath in 
% "Detection of signals by information theoretic criteria", IEEE
% Transactions on Acoustics, Speech and Signal Processing 33 (2) (1985)
% 387–392.

a   = zeros(min(p,n),1); 
MDL = zeros(min(p,n),1); 

for k=0:min(p,n)-1
    a(k+1) = 1/(p-k) * sum(lambda(k+1:p)); 
    MDL(k+1) = -(p-k)*n*(  1/(p-k) * sum(log(lambda(k+1:p))) - log(a(k+1)) ) ; 
    MDL(k+1) = MDL(k+1) + 1/2 * k * (2*p-k) * log(n);
end

[val, loc] = min(MDL); 

K = loc-1;

return; 
