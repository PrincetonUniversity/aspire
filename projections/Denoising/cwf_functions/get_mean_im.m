function [ mean_image] = get_mean_im( fn, mean_coeff, L0, R )
% Tejal Bhamre, Oct 2015
%INPUTS:
%fn: IFT of FB basis, stored in cell structure. cell number corresponds to the angular index.     
%mean_coeff: Coefficients of mean image 
%L0: Size of image
%R: Size of compact support
%OUTPUTS: mean image

L = 2*R;
tmp = fn{1};
tmp = reshape(tmp, L^2, size(tmp, 3));
mean_Im = reshape(tmp*mean_coeff, L, L);
mean_Im = real(mean_Im);

%Original image size
origin = floor(L0/2) + 1;
mean_image = zeros(L0);
mean_image(origin-R:origin+R-1, origin-R:origin+R-1) = mean_Im;

end

