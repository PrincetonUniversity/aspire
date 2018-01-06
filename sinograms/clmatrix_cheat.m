function [clmatrix,clcorr,cltheta]=clmatrix_cheat(rots,n_theta)
%
% Build common lines matrix using the true rotations corresponding to the
% projections orientations. Each projection has n_theta rays.
% clcorr is set to 1.0e-8.
%
% clmatrix contains the indices of the common lines between each pair of
% projections. cltheta contains the exact angle of the common lines and
% thus does not contain any discretization errors.
%
% Yoel Shkolnisky, October 2008.
%
% Revised, Y.S. June 2014.


N=size(rots,3);

clmatrix=zeros(N);   % common lines matrix
clcorr=zeros(N);     % correlation coefficient for ach common line
cltheta=zeros(N);    % angles of common line pairs

for k1=1:N-1
    R1=rots(:,:,k1);    
    for k2=k1+1:N               
        R2=rots(:,:,k2);
        [l_k1k2,l_k2k1]=commonline_R(R1,R2,n_theta);
        clmatrix(k1,k2)=l_k1k2+1;
        clmatrix(k2,k1)=l_k2k1+1;
        clcorr(k1,k2)=1.0e-8;
    end
end
