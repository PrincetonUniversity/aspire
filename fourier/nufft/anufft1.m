% ANUFFT1 Wrapper for adjoint non-uniform FFT (1D)
%
% Usage
%    sig = anufft1(sig_f, fourier_pts, sz);
%
% Input
%    sig_f: An Fourier transform calculated at the frequencies specified
%       by fourier_pts.
%    fourier_pts: The frequencies in Fourier space at which the adjoint Fourier
%       transform is to be calculated. These are in the form of a vector of
%       length K with values in the range [-pi, pi].
%
% Output
%    sig: The adjoint Fourier transform of sig_f at frequencies fourier_pts.

% See also
%    anudft1

function sig = anufft1(sig_f, fourier_pts, sz)
	persistent p_plan p_sz p_num_pts;

	epsilon = 1e-10;

	lib_code = pick_nufft_library(sz);

	num_pts = size(fourier_pts, 2);

	if lib_code == 3
		if ~isempty(p_plan) && all(p_sz==sz) && p_num_pts == num_pts
			plan = p_plan;
		else
			plan = nfft_init_1d(sz(1), num_pts);
		end

		nfft_set_x(plan, 1/(2*pi)*fourier_pts);
		nfft_precompute_psi(plan);
		nfft_set_f(plan, double(sig_f(:)));

		nfft_adjoint(plan);
		sig = nfft_get_f_hat(plan);

		if isempty(p_plan) || plan ~= p_plan
			if ~isempty(p_plan)
				nfft_finalize(p_plan);
			end
			p_plan = plan;
			p_sz = sz;
			p_num_pts = num_pts;
		end
	elseif lib_code == 2
		sig = num_pts*nufft1d1(num_pts, ...
			fourier_pts(1,:), ...
			double(sig_f(:)), 1, epsilon, sz(1));
	elseif lib_code == 1
        warning('NUFFT:directImplementation','Using direct (very slow) NUFFT. Call install_cims_nufft');
		sig = anudft1(sig_f, fourier_pts, sz);
	else
		error('invalid library code');
	end

	if isa(sig_f, 'single')
		sig = single(sig);
	end
end
