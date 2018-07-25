function [X,r,w,eigVecs,coef,sd_c] = matrixDenoise(Y,nv)
%MATRIXDENOISE: Perform asymptotically-optimal low rank singular-value shrinkage (according to Gavish-Donoho (2017))

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

% Y = Y/sqrt(n*nv);
[u,s,v] = svd(Y,'econ');
sd = diag(s);
% sd = sd(sd>(1+sqrt(beta)));
sd = sd(sd>(1+sqrt(beta))*sqrt(n*nv));
r = numel(sd);  
u = u(:,1:r);
v = v(:,1:r);
% sd_c = sqrt((sd.^2-beta-1).^2-4*beta)./sd;
sd_c = n*nv*sqrt((sd.^2/n/nv-beta-1).^2-4*beta)./sd;
w = sd_c./sd;

X = u*diag(sd_c)*v';
% X = sqrt(n*nv)*X;
if trans
    X = X.';
end

eigVecs = u;
coef = diag(sd)*v';

end

