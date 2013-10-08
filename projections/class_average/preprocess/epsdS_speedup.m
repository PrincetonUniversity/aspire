function [ P1d, P2 ] = epsdS_speedup( vol,samples_idx,biasflag, norm_flag, window_type,verbose )
%
% Estimate the 1D isotropic power spectrum function from image edges
%
%   vol - a tack of square pxp projections
%   samples_idx - is a list of pixel indices to use for autocorrelation estimation
%   max_d - the autocorrelation is estimated for distances up to max_d
%   biasflag - flag controls the bias of the estimate:
%    0    do not remove the mean from the samples before computing the
%         autocorrelation.
%    1    remove the mean
%    2    do not remove the mean but set the variance (autocorrelation at
%         distance 0) to be 1 (default).
%   R - autocorrelation function
%   norm_flag - normalize the esimate spectrum.
%   window_type - correlogram window type.
%  
%
% The power spectrum is esimated using correlogram method, mapping the correlation
%   over 2D circle of radius max_d and zero-padded to have an image of size
%   pxp.
%
%
%
% Other m-files required: epsdR.m, iwindow.m, cart2rad.m
% Subfunctions: none
% MAT-files required: none
%
% Author:  Yoel Shkolnisky, Mor Cohen
%________________________________________

% email: morcohe5@mail.tau.ac.il
% May 2011; Last revision: 04-May-2011
%
% See also: 


% Estimate the 1D isotropic autocorrelation function
[R,x]=epsdR_speedup(vol,samples_idx,biasflag,verbose);

if ~exist('window_type','var')
    window_type='boxcar'; % correlogram window default type
end

% Mapping the correlation over 2D circle of radius max_d and zero-padded 
% to have an image of size pxp
p=size(vol,1);
I = cart2rad(2*p+1);
[h, i,m] = unique(I);
R2 = zeros(length(h),1);
for n = 1:length(x) idx = h == x(n); R2(idx) = R(n); end

R2 = R2(m,:);
R2 = reshape(R2,size(I));
w = iwindow(2*p+1, window_type);

% FFT
P2 = cfft2(R2.*w);
P2=real(P2);
% Normalize the power spectrum
if norm_flag
    P2 = P2/norm(P2);
end   

% Demapping back to 1D function (since P2 is also "partially" isotropic - 
% assuming alising is negligible)
P1d = P2(i);

if (norm(imag(P1d))/real(P1d)) > eps
    error('error, psd suppose to be real function');
end 
    
P1d = real(P1d);

end

