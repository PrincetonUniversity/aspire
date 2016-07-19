% ANUFFT2 Wrapper for adjoint non-uniform FFT (2D)
%
% Usage
%    im = anufft2(im_f, fourier_pts, sz);
%
% Input
%    im_f: An array containing the Fourier transform of an image at certain
%       points.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 2-by-N array. These need to be in the
%       range [-pi, pi] in each dimension.
%    sz: The size of the resulting image.
%
% Output
%    image: The adjoint Fourier transform of im_f at the frequencies
%       fourier_pts.
%
% See also
%    anudft2

function im = anufft2(im_f, fourier_pts, sz)
	persistent p_plan p_sz p_num_pts;

	epsilon = 1e-10;

	lib_code = pick_nufft_library(sz);

	num_pts = size(fourier_pts, 2);

	if lib_code == 3
		if ~isempty(p_plan) && all(p_sz==sz) && p_num_pts == num_pts
			plan = p_plan;
		else
			plan = nfft_init_2d(sz(1), sz(2), num_pts);
		end

		nfft_set_x(plan, 1/(2*pi)*fourier_pts);
		nfft_precompute_psi(plan);
		nfft_set_f(plan, double(im_f(:)));

		nfft_adjoint(plan);
		im = nfft_get_f_hat(plan);
		im = permute(reshape(im, [sz(1) sz(2)]), [2 1]);

		if isempty(p_plan) || plan ~= p_plan
			if ~isempty(p_plan)
				nfft_finalize(p_plan);
			end
			p_plan = plan;
			p_sz = sz;
			p_num_pts = num_pts;
		end
	elseif lib_code == 2
		im = num_pts*nufft2d1(num_pts, ...
			fourier_pts(1,:), fourier_pts(2,:), ...
			double(im_f(:)), 1, epsilon, sz(1), sz(2));
		im = reshape(im, sz);
	elseif lib_code == 1
		im = anudft2(im_f, fourier_pts, sz);
	else
		error('invalid library code');
	end

	if isa(im_f, 'single')
		im = single(im);
	end
end
