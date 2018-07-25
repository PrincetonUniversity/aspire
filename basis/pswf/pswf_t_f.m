function [coeffs, PSWF_Nn_p] = pswf_t_f(images, L, beta, T, realFlag, PSWF_Nn_p)
% This is the forward PSWF transform - from images sampled on the Cartesian
% grid to the PSWF expansion coefficients.
%   Input:  images:     the images to be approximated, provided as a 3D array. The
%                       first two dimensions of the array correspond to a regular grid of (2L+1)x(2L+1)
%                       points, and the third dimension to the image index.
%           L:          Image resolution, defined such that there are 2L+1 samples
%                       for each dimension of the image.
%           beta:       Bandlimit ratio relative to the Nyquist rate, between 0 and 1. In general, the bandlimit is 
%                       c = beta*pi*L, therefore for beta = 1 there is no oversampling assumed. 
%                       This parameter controls the bandlimit of the PSWFs.
%           T:          Truncation parameter, between 0 and 1e6, which controls the length of the
%                       expansion and the approximation error. Smaller values (close to zero) guarantee smaller errors, yet longer expansions, and vice-versa. 
%                       Note: Due to numerical considerations, do not exceed 1e6.
%           realFlag:   Flag 0 or 1, which indiciates whether the images
%                       are assumed to be real valued. If yes, then only
%                       non-negative Fourier-mode indices of the PSWFs are used.
%           PSWF_Nn_p:  If empty, the function computes the
%                       PSWFs for the chosen parameters, uses them for the
%                       computation of coefficients, and returns them in PSWF_Nn_p. If not empty, the function will
%                       use the basis functions provided (use previousy pre-computed basis functions returned by this function). 
%                       Note that the basis functions are specific for parameters L, beta
%                       and T, and the function will not check whether the
%                       provided basis functtions correspond to the set of parameters.
%   Output: coeffs:     PSWF expansion coefficients. For real-valued images, only non-negative N indices are returned.
%           PSWF_Nn_p:  A struct which contains the PSWFs with their generation parameters, specifically:
%                           samples:    PSWF basis functions, samplned on the Cartesian grid inside the unit disk. Columns correspond to different indices. 
%                           Alpha_Nn:   PSWF associated eigenvalues.
%                           ang_freq:   Angular frequency ('N' indices) of the PSWFs.
%                           rad_freq:   Radial frequency ('n' indices) of the PSWFs.

%% Generate regular grid and basis functions (if needed)
x_1d_grid = -L:1:L;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= L);

if isempty(PSWF_Nn_p)
    x = x_2d_grid(points_inside_the_circle);
    y = y_2d_grid(points_inside_the_circle);
    [ PSWF_Nn_p.samples, PSWF_Nn_p.alpha_Nn, PSWF_Nn_p.ang_freq, PSWF_Nn_p.rad_freq ] = PSWF_gen_v3( L+1, L, beta, eps, T, x, y );
    PSWF_Nn_p.L = L;
    PSWF_Nn_p.beta = beta;
    PSWF_Nn_p.T = T;
end

%% compute expansion coefficients
nImages = size(images,3);
images_c = reshape(images,size(images,1)*size(images,2),nImages);    % Reshape from 3d array to a matrix
images_c = images_c(points_inside_the_circle(:),:);     % Take points inside the unit disk

coeffs = PSWF_Nn_p.samples' * images_c;
if (realFlag==0)
    coeffs = [PSWF_Nn_p.samples(:,PSWF_Nn_p.ang_freq>0).' * images_c ;coeffs];
end
    

end

