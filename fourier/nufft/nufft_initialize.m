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

	if lib_code == 3
		if numel(sz) == 1
			plan.nfft_plan_id = nfft_init_1d(sz(1), num_pts);
		elseif numel(sz) == 2
			plan.nfft_plan_id = nfft_init_2d(sz(1), sz(2), num_pts);
		elseif numel(sz) == 3
			plan.nfft_plan_id = nfft_init_3d(sz(1), sz(2), sz(3), num_pts);
		end
	elseif ~ismember(lib_code, [1 2 4])
		error('Invalid library code.');
	end
end
