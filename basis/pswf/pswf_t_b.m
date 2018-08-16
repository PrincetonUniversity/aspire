function images = pswf_t_b( coeffs, PSWF_Nn_p, realFlag )
% This is the backward PSWF transform - from PSWF expansion coefficients to images sampled on the Cartesian
% grid inside the unit disk.
%   Input:  coeffs:     PSWF expansion coefficients, in the format provided by the functions pswf_t_f.
%           PSWF_Nn_p:  PSWF struct, as provided by the functions pswf_t_f.
%           realFlag:   Flag 0 or 1, indiciating whether the images are assumed to be real valued.
%                       Note: Should be consistent with the provided coefficients

%% Compute image samples inside the unit disk
if (realFlag == 1)
    images_c = PSWF_Nn_p.samples(:,PSWF_Nn_p.ang_freq==0)*coeffs(PSWF_Nn_p.ang_freq==0,:) + 2*real(PSWF_Nn_p.samples(:,PSWF_Nn_p.ang_freq~=0)*coeffs(PSWF_Nn_p.ang_freq~=0,:));
else
    images_c = ([conj(PSWF_Nn_p.samples(:,PSWF_Nn_p.ang_freq>0)) PSWF_Nn_p.samples]*coeffs);
end

%% Save to 3D array
L = PSWF_Nn_p.L;
x_1d_grid = -L:1:L;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= L);

nImages = size(coeffs,2);
images = zeros(2*L+1,2*L+1,nImages);
for i=1:nImages
    tmp = zeros(2*L+1,2*L+1);
    tmp(points_inside_the_circle) = images_c(:,i);
    images(:,:,i) = tmp;
end

