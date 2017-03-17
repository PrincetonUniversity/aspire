% NUFFT1 Wrapper for non-uniform FFT (1D)
%
% Usage
%    sig_f = nufft1(sig, fourier_pts);
%
% Input
%    sig: A vector of length containing a signal.
%    fourier_pts: The frequencies in Fourier space at which the Fourier trans-
%       form is to be calculated. These are arranged as a vector of length K,
%       with values in the range [-pi, pi].
%
% Output
%    sig_f: The Fourier transform of sig at the frequencies fourier_pts.
%
% See also
%    nudft1

function sig_f = nufft1(sig, fourier_pts)
	persistent p_plan p_sz p_num_pts;

	epsilon = 1e-10;
	sz = size(sig);

	lib_code = pick_nufft_library(sz);

	num_pts = size(fourier_pts, 2);

	if lib_code == 3
		if ~isempty(p_plan) && all(p_sz==sz) && p_num_pts == num_pts
			plan = p_plan;
		else
			plan = nfft_init_1d(sz(1),num_pts);
		end

		nfft_set_x(plan, 1/(2*pi)*fourier_pts);
		nfft_precompute_psi(plan);
		nfft_set_f_hat(plan, sig);

		nfft_trafo(plan);
		sig_f = nfft_get_f(plan);

		if isempty(p_plan) || plan ~= p_plan
			if ~isempty(p_plan)
				nfft_finalize(p_plan);
			end
			p_plan = plan;
			p_sz = sz;
			p_num_pts = num_pts;
		end
	elseif lib_code == 2
		sig_f = nufft1d2(num_pts, ...
			fourier_pts(1,:), ...
			-1, epsilon, sz(1), double(sig(:)));
	elseif lib_code == 1
		sig_f = nudft1(sig, fourier_pts);
	else
		error('invalid library code');
	end

	if isa(sig, 'single')
		sig_f = single(sig_f);
	end
end
