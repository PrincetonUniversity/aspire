% ANUFFT3 Wrapper for adjoint non-uniform FFT (3D)
%
% Usage
%    vol = anufft3(vol_f, fourier_pts, sz);
%
% Input
%    vol_f: An array containing the Fourier transform of a volume at certain
%       points.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 3-by-N array. These need to be in the
%       range [-pi, pi] in each dimension.
%    sz: The size of the resulting volume.
%
% Output
%    vol: The adjoint Fourier transform of vol_f at the frequencies
%       fourier_pts.
%
% See also
%    anudft3

function vol = anufft3(vol_f, fourier_pts, sz)
	persistent p_plan p_sz p_num_pts;

	epsilon = 1e-10;

	lib_code = pick_nufft_library(sz);

	num_pts = size(fourier_pts, 2);

	if lib_code == 3
		if ~isempty(p_plan) && all(p_sz==sz) && p_num_pts == num_pts
			plan = p_plan;
		else
			plan = nfft_init_3d(sz(1), sz(2), sz(3), num_pts);
		end

		nfft_set_x(plan, 1/(2*pi)*fourier_pts);
		nfft_precompute_psi(plan);
		nfft_set_f(plan, double(vol_f(:)));

		nfft_adjoint(plan);
		vol = nfft_get_f_hat(plan);
		vol = permute(reshape(vol, [sz(1) sz(2) sz(3)]), [3 2 1]);

		if isempty(p_plan) || plan ~= p_plan
			if ~isempty(p_plan)
				nfft_finalize(p_plan);
			end
			p_plan = plan;
			p_sz = sz;
			p_num_pts = num_pts;
		end
	elseif lib_code == 2
		vol = num_pts*nufft3d1(num_pts, ...
			fourier_pts(1,:), fourier_pts(2,:), fourier_pts(3,:), ...
			double(vol_f(:)), 1, epsilon, sz(1), sz(2), sz(3));
		vol = reshape(vol, sz);
	elseif lib_code == 1
        warning('NUFFT:directImplementation','Using direct (very slow) NUFFT. Call install_cims_nufft');
		vol = anudft3(vol_f, fourier_pts, sz);
	else
		error('invalid library code');
	end

	if isa(vol_f, 'single')
		vol = single(vol);
	end
end
