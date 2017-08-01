% APPLY_MEAN_KERNEL Applies the mean kernel represented by convolution
%
% Usage
%    vol_basis = apply_mean_kernel(vol_basis, kernel_f, basis, mean_est_opt);
%
% Input
%    vol_basis: The volume to be convolved, stored in the basis coordinates.
%    kernel_f: The centered Fourier transform of the convolution kernel
%       representing the mean  projection-backprojection operator as obtained
%       from `cryo_mean_kernel_f`.
%    basis: A basis object corresponding to the basis used to store `vol_basis`.
%    mean_est_opt: An options structure. Currently no options are used.
%
% Output
%    vol_basis: The result of evaluating `vol_basis` in the given basis,
%       convolving with the kernel given by `kernel_f`, and backprojecting into
%       the basis.

function vol_basis = apply_mean_kernel(vol_basis, kernel_f, basis, mean_est_opt)
    vol = basis.evaluate(vol_basis);

    vol = cryo_conv_vol(vol, kernel_f);

    vol_basis = basis.evaluate_t(vol);
end
