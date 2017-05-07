function [pf,freqs]=cryo_pft(p,n_r,n_theta,precision)
%
% Compute the polar Fourier transform of projections with resolution n_r in
% the radial direction and resolution n_theta in the angular direction.
%
% If p is a volume, the function computes the polar Fourier transform of
% each slice in the volume seperately.
%
% Input parameters:
%   p          3D array. If p is a volume, then the last dimension
%              corresponds to the projection index: p(:,:,k) is projection
%              k.  Two first dimensions are x and y of the each projection.
%   n_r        Number of samples along each ray (in the radial direction).
%   n_theta    Angular resolution. Number of Fourier rays computed for each
%              projection.
%   precision  'single' or 'double'. Default is 'single'. The polar Fourier
%              samples for 'single' are computed to accuracy of 1.e-6.
%
% Output parameters:
%   pf      Polar Fourier transform of the input array. pf is a matrix of
%           with n_r rows and n_theta columns. Each column corresponds to
%           a fixed theta. The first column corresponds a theta=0. The
%           last column corresponds to theta nearest 2*pi. The first row
%           corresponds to r=0. The lase row correspond to r nearest pi.
%           If f is a volume with n slices, pf is a volume of size 
%           n_r x n_theta x n. The third index is the slice number; the
%           other two are as above.
%   freqs   Frequencies at which the polar Fourier samples were computed. A
%           matrix with n_rxn_theta rows and two columns (omega_x and
%           omega_y).
%          
% Revisions:
%   02/03/2009   Filename chaged from cryo_pft_v3.m to cryo_pft.m.
%   12/20/2009   omega0 changed from 2*pi/(2*n_r+1) to 2*pi/(2*n_r-1)
%
% Yoel Shkolnisky, January 2008.

if nargin<4
    precision='single';
end

n_proj=1;
if ndims(p)==3
    n_proj=size(p,3);
end
    
%n_uv=size(p,1);
omega0=2*pi/(2*n_r-1);
dtheta=2*pi/n_theta;

freqs=zeros(n_r*n_theta,2); % sampling points in the Fourier domain
for j=1:n_theta
    for k=1:n_r
        freqs((j-1)*n_r+k,:)=[(k-1)*omega0*sin((j-1)*dtheta),...
            (k-1)*omega0*cos((j-1)*dtheta)];
    end
end

%   freqs is the frequencies on [-pi,pi]x[-pi,pi] on which we sample the
%   Fourier transform of the projections. An array of size n_r*n_theta by 2
%   where each row corresponds to a frequnecy at which we sample the
%   Fourier transform of the projections. The first column is omega_x, the
%   second is omega_y. 


% precomputed interpolation weights once for the give polar grid. This is
% used below for computing the polar Fourier transform of all slices
precomp=nufft_t_2d_prepare(freqs,size(p,1),precision);

%t1=tic;
pf=zeros(n_r,n_theta,n_proj);
parfor k=1:n_proj
    tmp=p(:,:,k);
%    tmp=(2/n_uv)^2*nufft_t_v3(tmp,precomp);    
    tmp=nufft_t_2d_execute(tmp,precomp);    
    pf(:,:,k)=reshape(tmp,n_r,n_theta);   
end
%t1=toc(t1);

% The following code uses the nufft wrappers, which use either 'chemnitz'
% or'cims' (or defualt DFT). However, 'chemnitz' works only for even sizes,
% and for some of the tested size, 'cims' is 2 times slower than my code.
% So the following code is commented until the odd/even sizes would be
% fixed.
% t2=tic;
% pf2=zeros(n_r,n_theta,n_proj);
% for k=1:n_proj
%     tmp=p(:,:,k);    
%     tmp2=nufft2(tmp,-freqs.');
%     pf2(:,:,k)=reshape(tmp2,n_r,n_theta);   
% end
% t2=toc(t2);
% 
% err=norm(pf(:)-pf2(:))/norm(pf(:));
% fprintf('err = %e\n',err);
% fprintf('t1= %5.2f\n',t1);
% fprintf('t2= %5.2f\n',t2);
% fprintf('ratio = %5.2f\n',t1/t2);