% CIRCULARIZE_KERNEL_F Calculate circulant approximation of Fourier kernel
%
% Usage
%    kernel_circ_f = circularize_kernel_f(kernel_f);
%
% Input
%    kernel_f: The original kernel in Fourier of size 2N-by-2N-by-... that is
%       to be approximated by a circular convolution kernel of size
%       N-by-N-by-... .
%
% Output
%    kernel_circ_f: A convolution kernel in Fourier of size N-by-N-by-...
%       corresponding to the circular approximation of kernel_f that is
%       optimal in the Frobenius norm of the convolution operator.
%
% Note
%    This approximation corresponds to the "optimal preconditioner" described
%    by Tyrtyshnikov.

function kernel_circ_f = circularize_kernel_f(kernel_f)
    dims = ndims(kernel_f);

    kernel = mdim_fftshift(ifftn(kernel_f), 1:dims);

    kernel_circ = kernel;
    for dim = 1:dims
        kernel_circ = circularize_kernel_1d(kernel_circ, dim);
    end

    kernel_circ_f = fftn(mdim_ifftshift(kernel_circ, 1:dims));
end

function kernel_circ = circularize_kernel_1d(kernel, dim)
	sz = size(kernel);
	N = sz(dim)/2;

	% Set up subsref structure to extract "top" half.
	S.type = '()';
	S.subs = cell(1, length(sz));
	S.subs(:) = {':'};
	S.subs{dim} = 1:N;

	% Multiplier for weighted average.
	mult = reshape([0:N-1]/N, [ones(1,dim-1) N ones(1,length(sz)-dim-1)]);

	kernel_circ = bsxfun(@times, mult, subsref(kernel, S));

	% Do the "bottom" half.
	S.subs{dim} = N+1:2*N;

	mult = reshape([N:-1:1]/N, [ones(1,dim-1) N ones(1,length(sz)-dim-1)]);

	kernel_circ = kernel_circ + bsxfun(@times, mult, subsref(kernel, S));

	% Shift to make sure zero frequency is in the center.
	kernel_circ = fftshift(kernel_circ, dim);
end
