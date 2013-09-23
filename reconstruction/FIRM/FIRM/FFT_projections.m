function [ projs_fourier ] = FFT_projections( projections, shifts )
% FFT on projections, and add phase to the Fourier slices according to the
% shifts
% Input: 
%
%   projections: stack of projections, of size n x n x n_proj, where
%   n is the size of a projection, and n_proj is the number of the
%   projections.
%   
%   shifts: the translation parameters for the projections. Default: [].
%
% Output:
%   projs_fourier: the square Fourier slices, a stack of size n x n x n_proj.
%
% Lanhui Wang Feb 10, 2012

%% FFT on the projections
n_proj=size(projections,3);
n=size(projections,1);
projs_fourier=fftshift(fft2(ifftshift(projections)));

%% Correct the Fourier coefficients according to the shift information
if nargin<2
    shifts=[];
    fprintf('No shift correction will be performed.\n');
else
    if mod(n,2)==0
        [omega_x,omega_y]=ndgrid(-n/2:(n/2-1),-n/2:(n/2-1));
    else
        [omega_x,omega_y]=ndgrid(-(n-1)/2:(n-1)/2,-(n-1)/2:(n-1)/2);
    end
    omega_x=-2*pi.*omega_x/n; omega_y=-2*pi.*omega_y/n;
    for k=1:n_proj
        
        p_fourier=projs_fourier(:,:,k);
        phase_x=omega_x.*shifts(k,1);
        phase_y=omega_y.*shifts(k,2);
        p_fourier=p_fourier.*exp(sqrt(-1)*(phase_x+phase_y));
        projs_fourier(:,:,k)=p_fourier;
        
    end
end

end