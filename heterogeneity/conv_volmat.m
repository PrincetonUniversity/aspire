% CONV_VOLMAT Convolve volume matrix with kernel
%
% Usage
%    volmat = conv_volmat(volmat, kernel_f);
%
% Input
%    volmat: An N-by-...-by-N (6 dimensions) volume matrix to be convolved.
%    kernel_f: The Fourier transform of the cubic matrix convolution kernel.
%        Must be larger than volmat. This must be a non-centered Fourier
%        transform.
%
% Output
%    volmat: The original volume matrix convolved by the kernel with the same
%       dimensions as before.

function volmat = conv_volmat(volmat, kernel_f)
    sz_volmat = size(volmat);
    N = sz_volmat(1);

    if any(sz_volmat(1:6)~=N)
        error('volume matrix must be cubic and square');
    end

    % TODO: Deal with rolled dimensions

    is_singleton = numel(sz_volmat)==6;

    sz_ker = size(kernel_f);
    N_ker = sz_ker(1);

    if any(sz_ker~=N_ker)
        error('convolution kernel must be cubic and square');
    end

    if ~isoctave
        % NOTE: Order is important here. It's about 20% faster to run from
        % 1 through 6 compared with 6 through 1.
        volmat = fft(volmat, N_ker, 1);
        volmat = fft(volmat, N_ker, 2);
        volmat = fft(volmat, N_ker, 3);
        volmat = fft(volmat, N_ker, 4);
        volmat = fft(volmat, N_ker, 5);
        volmat = fft(volmat, N_ker, 6);

        volmat = bsxfun(@times, volmat, kernel_f);

        % NOTE: Again, the order here is important. Also, we can't use
        % 'symmetric' option here since the partial IFFTs are not real.
        volmat = ifft(volmat, [], 6);
        volmat = volmat(:,:,:,:,:,1:N,:);
        volmat = ifft(volmat, [], 5);
        volmat = volmat(:,:,:,:,1:N,:,:);
        volmat = ifft(volmat, [], 4);
        volmat = volmat(:,:,:,1:N,:,:,:);
        volmat = ifft(volmat, [], 3);
        volmat = volmat(:,:,1:N,:,:,:,:);
        volmat = ifft(volmat, [], 2);
        volmat = volmat(:,1:N,:,:,:,:,:);
        volmat = ifft(volmat, [], 1);
        volmat = volmat(1:N,:,:,:,:,:,:);
    else
        if is_singleton
            volmat = fftn(volmat, N_ker*ones(1, 6));
        else
            volmat = fft(volmat, N_ker, 1);
            volmat = fft(volmat, N_ker, 2);
            volmat = fft(volmat, N_ker, 3);
            volmat = fft(volmat, N_ker, 4);
            volmat = fft(volmat, N_ker, 5);
            volmat = fft(volmat, N_ker, 6);
        end

        volmat = bsxfun(@times, volmat, kernel_f);

        if is_singleton
            volmat = ifftn(volmat);
            volmat = volmat(1:N,1:N,1:N,1:N,1:N,1:N);
        else
            volmat = ifft(volmat, [], 6);
            volmat = volmat(:,:,:,:,:,1:N,:);
            volmat = ifft(volmat, [], 5);
            volmat = volmat(:,:,:,:,1:N,:,:);
            volmat = ifft(volmat, [], 4);
            volmat = volmat(:,:,:,1:N,:,:,:);
            volmat = ifft(volmat, [], 3);
            volmat = volmat(:,:,1:N,:,:,:,:);
            volmat = ifft(volmat, [], 2);
            volmat = volmat(:,1:N,:,:,:,:,:);
            volmat = ifft(volmat, [], 1);
            volmat = volmat(1:N,:,:,:,:,:,:);
        end
    end
end

