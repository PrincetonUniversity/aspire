function [X,r,w,eigVecs,coef,sd_c,sd] = matrixDenoise_v2(Y,nv,useReflections)
%MATRIXDENOISE: Perform asymptotically-optimal low-rank matrix denoising by singular value shrinkage (according to Gavish-Donoho (2017))

[sizeX,sizeY] = size(Y);
if (sizeX>sizeY)    % Define matrix to be wide
    Y = Y.';
    m = sizeY; 
    n = sizeX;
    trans = true;
else
    m = sizeX;
    n = sizeY;
    trans = false;
end
beta = m/n;

if ~useReflections      % Don't include reflections -- Compute standard SVD
    [u,s,v] = svd(Y,'econ');
    sd = diag(s);
else                    % Include reflections -- Obtain singular components from the real part of the covariance
    [u,d] = eig(real(Y*Y')/n);
    [d,idx] =sort(diag(d),'descend'); u=u(:,idx);
    sd = sqrt(d*n);
    v = u'*Y; 
    v = v./(sum(abs(v).^2,2)*ones(1,n));
    v = v';
end

sd = sd(sd>(1+sqrt(beta))*sqrt(n*nv));  % Baik-Ben Arous-Peche transition
r = numel(sd);  
u = u(:,1:r);
v = v(:,1:r);
sd_c = n*nv*sqrt((sd.^2/n/nv-beta-1).^2-4*beta)./sd;    % Shrunk singular values
w = sd_c./sd;   % Weight = ratio between shrunk and original singular values

X = u*diag(sd_c)*v';    % Denoised matrix

if trans
    X = X.';
end

eigVecs = u;        % Eigenvectors for steerable PCA
coef = diag(sd)*v'; % Coefficients for steerable PCA (no shrinkage)

end

