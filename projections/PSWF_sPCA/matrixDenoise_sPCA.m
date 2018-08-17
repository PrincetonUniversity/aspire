function [X,r,w,eigVecs,coef,coef_den,lambda] = matrixDenoise_sPCA(Y,nv,useReflections)
%MATRIXDENOISE: Perform asymptotically-optimal low-rank matrix denoising by singular value shrinkage (according to Gavish-Donoho (2017))

[m,n] = size(Y);
if ~useReflections
    [u,s,~] = svd(Y,'econ');
    s = diag(s);
else        % Use reflections -- Implicitly augment dataset with reflections
    n=2*n;
    [u,d] = eig(2*real(Y*Y')/n);
    [d,idx] =sort(diag(d),'descend'); u=u(:,idx);
    s = sqrt(d*n);
end

beta = m/n;

s = s(s>(1+sqrt(beta))*sqrt(n*nv));  % Baik-Ben Arous-Peche transition
r = numel(s);  
u = u(:,1:r);
s_c = n*nv*sqrt((s.^2/n/nv-beta-1).^2-4*beta)./s;    % Shrunk singular values
w = s_c./s;       % Weight = ratio between shrunk and original singular values

lambda = s.^2/n;    % Eigenvalues for steerable PCA
eigVecs = u;        % Eigenvectors for steerable PCA
coef = u'*Y;        % Coefficients for steerable PCA (no shrinkage)
% coef_den = bsxfun(@times,w,coef);
coef_den = diag(w)*coef;
X = u*coef_den; % Denoised matrix 

end

