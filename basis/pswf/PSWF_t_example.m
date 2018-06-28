% clc;
close all;

%% Load example dataset 
load 70S_proj_10_s4     % Example with 10 projection images

%% Define parameters
L = 64;         % Image resolution: (2L+1)x(2L+1) samples on a Cartesian grid
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values stand for greater oversampling
c = beta*pi*L;  % Bandlimit 
T = 1e-3;       % Truncation parameter
realFlag = 1;   % Flag indicating whether the data is assumed to be real-valued

% For all parameters - see more details in the function pswf_t_f. Also, see the
% paper "Approximation scheme for essentially bandlimited and space-concentrated functions on a disk" by Boris Landa and Yoel Shkolnisky.

%% Obtain expansion coefficients
[coeffs, PSWF_Nn_p] = pswf_t_f(projections, L, beta, T, realFlag, []);          % Obtain coefficients while evaluating PSWFs
% [coeffs, PSWF_Nn_p] = pswf_t_f(projections, L, beta, T, realFlag, PSWF_Nn_p); % Obtain coefficients using precomputed PSWFs (first run this function with an empty array instead of PSWF_Nn_p)

%% Reconstruct images from coefficients
reconstructedImg = pswf_t_b( coeffs, PSWF_Nn_p, realFlag );

%% Compute errors and display reconstructed images
x_1d_grid = -L:1:L;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= L);
points_outside_the_circle = (r_2d_grid >= L);

err = reconstructedImg - projections;
e = mean(abs(err).^2,3);
e = sum(e(points_inside_the_circle));
p = mean(abs(projections).^2,3);
p = sum(p(points_inside_the_circle));
display(['Ratio between avarege squared error and image power: ', num2str(e/p)])

figure; subplot(1,3,1); imagesc(projections(:,:,1)); colormap(gray); title('Original image') 
subplot(1,3,2); imagesc(reconstructedImg(:,:,1)); colormap(gray); title('Approximated image')
err2Disp = abs(err(:,:,1)); err2Disp(points_outside_the_circle) = 0; % We can ignore values outside the disk
subplot(1,3,3); imagesc(err2Disp); colormap(gray); title('Error pattern')
