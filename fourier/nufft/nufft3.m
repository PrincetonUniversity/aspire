% NUFFT3 Wrapper for non-uniform FFT (3D)
%
% Usage
%    vol_f = nufft3(vol, fourier_pts);
%
% Input
%    vol: An N-by-N-by-N array containing the voxel structure of a volume.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as a 3-by-N array. These need to be in the
%       range [-pi, pi] in each dimension.
%
% Output
%    vol_f: The Fourier transform of vol at the frequencies fourier_pts.
%
% See also
%    nudft3

function vol_f = nufft3(vol, fourier_pts)
	persistent p_plan p_sz p_num_pts;

	epsilon = 1e-10;
	sz = size(vol);

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
		vol = reshape(permute(double(vol), [3 2 1]), prod(sz), 1);
		nfft_set_f_hat(plan, vol);

		nfft_trafo(plan);
		vol_f = nfft_get_f(plan);

		if isempty(p_plan) || plan ~= p_plan
			if ~isempty(p_plan)
				nfft_finalize(p_plan);
			end
			p_plan = plan;
			p_sz = sz;
			p_num_pts = num_pts;
		end
	elseif lib_code == 2
		vol_f = nufft3d2(num_pts, ...
			fourier_pts(1,:), fourier_pts(2,:), fourier_pts(3,:), ...
			-1, epsilon, sz(1), sz(2), sz(3), double(vol(:)));
	elseif lib_code == 1
		vol_f = nudft3(vol, fourier_pts);
	else
		error('invalid library code');
	end

	if isa(vol, 'single')
		vol_f = single(vol_f);
	end
end
