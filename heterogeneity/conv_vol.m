% CONV_VOL Convolve volume with kernel
%
% Usage
%    vol = conv_vol(vol, kernel_f);
%
% Input
%    vol: An N-by-N-by-N volume to be convolved.
%    kernel_f: The Fourier transform of the cubic convolution kernel. Must be
%       larger than vol. This must be a non-centered Fourier transform.
%
% Output
%    vol: The original volume convolved by the kernel with the same dimensions
%       as before.

function x = conv_vol(x, kernel_mu_f)
    N = size(x, 1);

    if any([size(x, 1) size(x, 2) size(x, 3)]~=N)
        error('volumes must be cubic');
    end

    N_ker = size(kernel_mu_f, 1);

    if any(size(kernel_mu_f)~=N_ker)
        error('convolution kernel must be cubic');
    end

    if ~isoctave
        padded_fft3 = @(x)(fft(fft(fft(x, N_ker, 3), ...
             N_ker, 2), N_ker, 1));

        x = ifftn(bsxfun(@times, padded_fft3(x), kernel_mu_f), [], 'symmetric');
    else
        padded_fft3 = @(x)(fftn(x, N_ker*ones(1, 3)));

        x = ifftn(bsxfun(@times, padded_fft3(x), kernel_mu_f));
    end
    x = x(1:N,1:N,1:N,:);
end
