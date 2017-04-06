function [shifts, corr, averages, norm_variance] = align_main_inmem( ...
    images, angle, class_VDM, refl, FBsPCA_data, k, max_shifts, list_recon, recon_sPCA)

% Function for aligning images with its k nearest neighbors to generate
% class averages. Wrapper for the align_main function providing in-memory
% functionality.
%
%   Usage:
%       [shifts, corr, averages, norm_variance] = align_main_inmem(images, ...
%          angle, class_VDM, refl, FBsPCA_data, k, max_shifts, list_recon);
%   Input:
%       images: LxLxP matrix. P projection images of size LxL pixels.
%       angle: Pxl (l>=k) matrix. Rotational alignment angle
%       class_VDM: Pxl matrix. Nearest neighbor list
%       refl: Pxl matrix. 1: no reflection. 2: reflection
%       FBsPCA_data: Fourier-Bessel steerable PCA data with r_max, UU,
%       Coeff, Mean, Freqs
%       k: number of nearest neighbors for class averages
%       max_shifts: maximum number of pixels to check for shift
%       list_recon: indices for images to compute class averages
%       recon_sPCA: reconstructed images after sPCA step
%   Output:
%       shifts: Pxk matrix. Relative shifts for k nearest neighbors
%       corr: Pxk matrix. Normalized cross correlation of each image with
%       its k nearest neighbors
%       average: LxLxP matrix. Class averages
%       norm_variance: compute the variance of each class averages.
%
% Joakim Andén Nov 2016

    data.dim = size(images);
    data.getImage = @(s)(images(:,:,s));

    tmpdir = tempname;
    mkdir(tmpdir);

    [shifts, corr, ave_filename, norm_variance] = align_main(data, angle, ...
        class_VDM, refl, FBsPCA_data, k, max_shifts, list_recon, recon_sPCA, tmpdir);

    averages = ReadMRC(ave_filename);

    delete(ave_filename);
    if isoctave()
        confirm_recursive_rmdir(false, 'local');
    end
    rmdir(tmpdir, 's');
end
