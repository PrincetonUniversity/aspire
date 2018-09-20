% NUFFT_INITIALIZE Initialize NUFFT plan
%
% Usage
%    plan = nufft_initialize(sz, num_pts, nufft_opt);
%
% Input
%    sz: A size vector of length 1, 2, or 3.
%    num_pts: The number of Fourier nodes.
%    nufft_opt: A struct containing the fields:
%       - epsilon: The desired precision of the NUFFT (default 1e-15).
%
% Output
%    plan: An NUFFT plan.

function plan = nufft_initialize(sz, num_pts, nufft_opt)
	if nargin < 3
		nufft_opt = [];
	end

	if numel(sz) < 1 || numel(sz) > 3
		error('Input ''sz'' vector must be of length 1, 2, or 3.');
	end

	if num_pts < 1 || ~isfinite(num_pts) || floor(num_pts) ~= num_pts
		error('Input ''num_pts'' must be a positive finite integer.')
	end

	nufft_opt = fill_struct(nufft_opt, ...
		'epsilon', 1e-15, ...
		'num_threads', 0);

	lib_code = pick_nufft_library(sz);

	plan = struct();

	plan.lib_code = lib_code;
	plan.sz = sz;
	plan.num_pts = num_pts;

	plan.epsilon = nufft_opt.epsilon;
	plan.num_threads = nufft_opt.num_threads;

	if lib_code == 3
		m = epsilon_to_nfft_cutoff(plan.epsilon);

		NFFT_SORT_NODES = bitshift(uint32(1), 11);
		FFTW_DESTROY_INPUT = bitshift(uint32(1), 0);

		nfft_flags = uint32(0);
		nfft_flags = bitor(nfft_flags, uint32(PRE_PHI_HUT));
		nfft_flags = bitor(nfft_flags, uint32(PRE_PSI));
		nfft_flags = bitor(nfft_flags, uint32(FFT_OUT_OF_PLACE));

		if numel(sz) > 1
			nfft_flags = bitor(nfft_flags, uint32(NFFT_SORT_NODES));
			nfft_flags = ...
				bitor(nfft_flags, uint32(NFFT_OMP_BLOCKWISE_ADJOINT));
		end

		fftw_flags = uint32(0);
		fftw_flags = bitor(fftw_flags, uint32(FFTW_ESTIMATE));
		fftw_flags = bitor(fftw_flags, uint32(FFTW_DESTROY_INPUT));

		sz1 = 2*pow2(nextpow2(sz));

		if numel(sz) == 1
			plan.nfft_plan_id = nfft_init_guru(1, sz(1), ...
				num_pts, sz1(1), m, nfft_flags, fftw_flags);
		elseif numel(sz) == 2
			plan.nfft_plan_id = nfft_init_guru(2, sz(1), sz(2), ...
				num_pts, sz1(1), sz1(2), m, nfft_flags, fftw_flags);
		elseif numel(sz) == 3
			plan.nfft_plan_id = nfft_init_guru(3, sz(1), sz(2), sz(3), ...
				num_pts, sz1(1), sz1(2), sz1(3), m, nfft_flags, fftw_flags);
		end
	elseif ~ismember(lib_code, [1 2 4 5])
		error('Invalid library code.');
	end
end

function m = epsilon_to_nfft_cutoff(epsilon)
	% NOTE: These are obtained empirically. Should have a theoretical
	% derivation.
	rel_errs = [6e-2 2e-3 2e-5 2e-7 3e-9 4e-11 4e-13 0];

	m = find(epsilon > rel_errs, 1);
end
