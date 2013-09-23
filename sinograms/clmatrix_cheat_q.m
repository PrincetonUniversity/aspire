function [clmatrix,clcorr]=clmatrix_cheat_q(q,n_theta)
%
% Build common lines matrix using the true quaternions corresponding to the
% projections orientations. Each projection has n_theta rays.
% clcorr is set to 1.0e-8.
%
% Yoel Shkolnisky, October 2008.

N=size(q,2);

clmatrix=zeros(N);   % common lines matrix
clcorr=zeros(N);     % correlation coefficient for ach common line

for k1=1:N-1
    for k2=k1+1:N
        [idx1,idx2]=commonline_q(q,k1,k2,n_theta);
        clmatrix(k1,k2)=idx1+1;
        clmatrix(k2,k1)=idx2+1;
        clcorr(k1,k2)=1.0e-8;
    end
end