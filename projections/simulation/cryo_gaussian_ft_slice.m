function [y,nf]=cryo_gaussian_ft_slice(def,n,rmax,q)
%
% Sample the analytic Fourier transform of the Gaussian phantom.
%
% Input parameteres:
%     def   Filename of parameters file of the phantom.
%     n     Dimension of each projection.
%     rmax  Physical dimensions of each projection. Each Fourier
%     transformed projection has n samples in each dimension between -rmax
%     and rmax. 
%     q     Array of quaternions.
%
% Output parameters:
%     y     3D array of size n x n x length(q). y(:,:,k) contains the
%           samples of the analytic the Fourier transform of the projection
%           whose direction determined by the k'th quaternion in q.
%     nf    Normalization factor. The icfft2 of the returned Fourier
%           transforms should be multiplied by nf to convert them into
%           projection images. See test_cryo_project for an example.
%
% Example:
%     voldef='C1_params';
%     q=qrand(1);
%     rmax=2; % I used larger rmax just to show that other rmax work too.
%     n=257; % Should be large enough for the given rmax.
%     [p1hat,nf]=cryo_gaussian_ft(voldef,q,n,rmax);
%     p1=icfft2(p1hat).*nf; 
%     assert(max(abs(imag(p1(:))))<1.0e-14)
%     p1=real(p1);
% 
%     p2=cryo_project_gaussian(n,rmax,q,voldef);
%     subplot(1,3,1); imagesc(p1); colorbar; axis image;
%     subplot(1,3,2); imagesc(p2); colorbar; axis image;
%     subplot(1,3,3); imagesc(p2-p1); colorbar; axis image;
%     err=norm(p1(:)-p2(:))/norm(p1(:));
%     disp(err); % Error should be close to machine precision.
%
% Yoel Shkolnisky, August 2013.

if mod(n,2)==1
    rng=-(n-1)/2:(n-1)/2;
    t=rng./ ((n-1)/2)*rmax;
else
    rng=-n/2+1/2:n/2-1/2;
    t=rng./(n/2)*rmax;
end


[I,J]=ndgrid(rng,rng); 
I=I(:); J=J(:);
T0=(t(2)-t(1))*n;
omega0=2*pi/T0;

nproj=size(q,2);
phat=zeros(n,n,nproj);

for k=1:nproj
    R=q_to_rot(q);
    Rt=R.';
    n_x= Rt(:,1);
    n_y= Rt(:,2);
    P = I * n_x' + J * n_y';    
    omega= omega0.*P;
    
    y=cryo_gaussian_ft(def,omega);
    y=reshape(y,n,n);
    phat(:,:,k)=y;
end
nf=(n/T0)^2; % For n odd, this is actually equal to ((n-1)/(2*rmax))^2, 
             % from which we can get a simple expression for T0 above. 
             % For n even, the constant is (n/(2*rmax))^2.

