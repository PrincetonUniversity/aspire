% NUFFT2 Wrapper for non-uniform FFT (2D)
%
% Usage
%    im_f = nufft2(im, fourier_pts);
%
% Input
%    im: An N-by-N array containing the pixel structure of an image.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 2-by-K array. These need to be in the
%       range [-pi, pi] in each dimension.
%
% Output
%    im_f: The Fourier transform of im at the frequencies fourier_pts.
%
% See also
%    nudft2

function im_f = nufft2(im, fourier_pts)
	p = nufft_initialize(size(im), size(fourier_pts, 2));

<<<<<<< HEAD
	epsilon = 1e-8;
	sz = size(im);
=======
	p = nufft_set_points(p, fourier_pts);
>>>>>>> master

	im_f = nufft_transform(p, im);

<<<<<<< HEAD
	num_pts = size(fourier_pts, 2);

	if lib_code == 3
		if ~isempty(p_plan) && all(p_sz==sz) && p_num_pts == num_pts
			plan = p_plan;
		else
			plan = nfft_init_2d(sz(1), sz(2), num_pts);
		end

		nfft_set_x(plan, 1/(2*pi)*fourier_pts);
		nfft_precompute_psi(plan);
		im = reshape(permute(double(im), [2 1]), prod(sz), 1);
		nfft_set_f_hat(plan, im);

		nfft_trafo(plan);
		im_f = nfft_get_f(plan);

		if isempty(p_plan) || plan ~= p_plan
			if ~isempty(p_plan)
				nfft_finalize(p_plan);
			end
			p_plan = plan;
			p_sz = sz;
			p_num_pts = num_pts;
		end
	elseif lib_code == 2
		im_f = nufft2d2(num_pts, ...
			fourier_pts(1,:), fourier_pts(2,:), ...
			-1, epsilon, sz(1), sz(2), double(im(:)));
	elseif lib_code == 1
        warning('NUFFT:directImplementation','Using direct (very slow) NUFFT. Call install_cims_nufft');
		im_f = nudft2(im, fourier_pts);
	else
		error('invalid library code');
	end

	if isa(im, 'single')
		im_f = single(im_f);
	end
=======
	nufft_finalize(p);
>>>>>>> master
end
