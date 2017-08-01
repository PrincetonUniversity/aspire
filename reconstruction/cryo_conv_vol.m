% CRYO_CONV_VOL Convolve volume with kernel
%
% Usage
%    vol = cryo_conv_vol(vol, kernel_f);
%
% Input
%    vol: An N-by-N-by-N-by-... array of volumes to be convolved.
%    kernel_f: The Fourier transform of the cubic convolution kernel. Must be
%       larger than vol in the first three dimensions. This must also be a
%       centered Fourier transform.
%
% Output
%    vol: The original volumes convolved by the kernel with the same dimensions
%       as before.

function x = cryo_conv_vol(x, kernel_f)
    N = size(x, 1);

    [x, sz_roll] = unroll_dim(x, 4);

    if any([size(x, 1) size(x, 2) size(x, 3)] ~= N)
        error('Volumes in `x` must be cubic.');
    end

    is_singleton = (numel(size(x)) == 3);

    N_ker = size(kernel_f, 1);

    if any(size(kernel_f) ~= N_ker)
        error('Convolution kernel `kernel_f` must be cubic.');
    end

    kernel_f = ifftshift(ifftshift(ifftshift(kernel_f, 1), 2), 3);

    if is_singleton
        x = fftn(x, N_ker*ones(1, 3));
    else
        x = fft(x, N_ker, 1);
        x = fft(x, N_ker, 2);
        x = fft(x, N_ker, 3);
    end

    x = bsxfun(@times, x, kernel_f);

    if is_singleton
        if ~isoctave
            x = ifftn(x, [], 'symmetric');
        else
            x = ifftn(x);
        end 

        x = x(1:N,1:N,1:N,:);
    else
        x = ifft(x, [], 1);
        x = x(1:N,:,:,:);
        x = ifft(x, [], 2);
        x = x(:,1:N,:,:);
        x = ifft(x, [], 3);
        x = x(:,:,1:N,:);
    end

    x = real(x);

    x = roll_dim(x, sz_roll);
end
