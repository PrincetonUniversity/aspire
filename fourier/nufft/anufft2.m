% ANUFFT2 Wrapper for adjoint non-uniform FFT (2D)
%
% Usage
%    im = anufft2(im_f, fourier_pts, sz);
%
% Input
%    im_f: An array containing the Fourier transform of an image at K
%       frequencies.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 2-by-K array. These need to be in the
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
	if numel(sz) ~= 2
		error('Input ''sz'' must have two elements.');
	end

	p = nufft_initialize(sz, size(fourier_pts, 2));

	p = nufft_set_points(p, fourier_pts);

<<<<<<< HEAD
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
        warning('NUFFT:directImplementation','Using direct (very slow) NUFFT. Call install_cims_nufft');
		im = anudft2(im_f, fourier_pts, sz);
	else
		error('invalid library code');
	end
=======
	im = nufft_adjoint(p, im_f);
>>>>>>> master

	nufft_finalize(p);
end
