% NUFFT_INITIALIZE Initialize NUFFT plan
%
% Usage
%    plan = nufft_initialize(sz, num_pts);
%
% Input
%    sz: A size vector of length 1, 2, or 3.
%    num_pts: The number of Fourier nodes.
%
% Output
%    plan: An NUFFT plan.

function plan = nufft_initialize(sz, num_pts)
	if numel(sz) < 1 || numel(sz) > 3
		error('Input ''sz'' vector must be of length 1, 2, or 3.');
	end

	if num_pts < 1 || ~isfinite(num_pts) || floor(num_pts) ~= num_pts
		error('Input ''num_pts'' must be a positive finite integer.')
	end

	lib_code = pick_nufft_library(sz);

	plan = struct();

	plan.lib_code = lib_code;
	plan.sz = sz;
	plan.num_pts = num_pts;

	plan.epsilon = 1e-15;

	if lib_code == 3
		m = epsilon_to_nfft_cutoff(plan.epsilon);

		nfft_flags = PRE_PHI_HUT | PRE_PSI | FFT_OUT_OF_PLACE;
		fftw_flags = FFTW_ESTIMATE;

		if numel(sz) == 1
			plan.nfft_plan_id = nfft_init_guru(1, sz(1), ...
				num_pts, 2*sz(1), m, nfft_flags, fftw_flags);
		elseif numel(sz) == 2
			plan.nfft_plan_id = nfft_init_guru(2, sz(1), sz(2), ...
				num_pts, 2*sz(1), 2*sz(2), m, nfft_flags, fftw_flags);
		elseif numel(sz) == 3
			plan.nfft_plan_id = nfft_init_guru(3, sz(1), sz(2), sz(3), ...
				num_pts, 2*sz(1), 2*sz(2), 2*sz(3), m, nfft_flags, fftw_flags);
		end
	elseif ~ismember(lib_code, [1 2 4])
		error('Invalid library code.');
	end
end

function m = epsilon_to_nfft_cutoff(epsilon)
	% NOTE: These are obtained empirically. Should have a theoretical
	% derivation.
	rel_errs = [6e-2 2e-3 2e-5 2e-7 3e-9 4e-11 4e-13 0];

	m = find(epsilon > rel_errs, 1);
end
