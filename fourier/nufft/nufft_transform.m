% NUFFT_TRANSFORM Apply transform according to NUFFT plan
%
% Usage
%    sig_f = nufft_transform(plan, sig);
%
% Input
%    plan: An NUFFT plan.
%    sig: A signal (1D, 2D, or 3D).
%
% Output
%    sig_f: Non-uniform Fourier transform of sig.

function sig_f = nufft_transform(plan, sig)
	if ~isfield(plan, 'lib_code') || ~isfield(plan, 'sz') || ...
		~isfield(plan, 'num_pts')
		error('Input ''plan'' is not a valid NUFFT plan.');
	end

	if ~isfield(plan, 'fourier_pts')
		error('Plan has not been initialized with Fourier points.');
	end

	dims = numel(plan.sz);

	sig_sz = size(sig);

	if dims == 1 && sig_sz(2) == 1
		sig_sz = sig_sz(1);
	end

	if numel(sig_sz) ~= dims || any(sig_sz ~= plan.sz)
		error('Input ''sig'' must be of size plan.sz.');
	end

	precision = class(sig);

	if plan.lib_code == 1
		if dims == 1
			sig_f = nudft1(sig, plan.fourier_pts);
		elseif dims == 2
			sig_f = nudft2(sig, plan.fourier_pts);
		elseif dims == 3
			sig_f = nudft3(sig, plan.fourier_pts);
		end
	elseif plan.lib_code == 2
		if strcmp(precision, 'double')
			epsilon = 1e-16;
		elseif strcmp(precision, 'single')
			epsilon = 1e-8;
		end

		sig = double(sig(:));

		if dims == 1
			sig_f = nufft1d2(plan.num_pts, ...
				plan.fourier_pts(1,:), ...
				-1, epsilon, plan.sz(1), sig);
		elseif dims == 2
			sig_f = nufft2d2(plan.num_pts, ...
				plan.fourier_pts(1,:), ...
				plan.fourier_pts(2,:), ...
				-1, epsilon, ...
				plan.sz(1), plan.sz(2), sig);
		elseif dims == 3
			sig_f = nufft3d2(plan.num_pts, ...
				plan.fourier_pts(1,:), ...
				plan.fourier_pts(2,:), ...
				plan.fourier_pts(3,:), ...
				-1, epsilon, ...
				plan.sz(1), plan.sz(2), plan.sz(3), sig);
		end
	elseif plan.lib_code == 3
		if ~isfield(plan, 'nfft_plan_id')
			error('Input ''plan'' is not a valid NUFFT plan.');
		end

		sig = double(sig);

		if dims == 2
			sig = reshape(permute(sig, [2 1]), prod(plan.sz), 1);
		elseif dims == 3
			sig = reshape(permute(sig, [3 2 1]), prod(plan.sz), 1);
		end

		nfft_set_f_hat(plan.nfft_plan_id, sig);
		nfft_trafo(plan.nfft_plan_id);
		sig_f = nfft_get_f(plan.nfft_plan_id);
	elseif plan.lib_code == 4
		if strcmp(precision, 'double')
			epsilon = 1e-16;
		elseif strcmp(precision, 'single')
			epsilon = 1e-8;
		end

		if dims == 1
			sig_f = finufft1d2( ...
				plan.fourier_pts(1,:), ...
				-1, epsilon, sig);
		elseif dims == 2
			sig_f = finufft2d2( ...
				plan.fourier_pts(1,:), ...
				plan.fourier_pts(2,:), ...
				-1, epsilon, sig);
		elseif dims == 3
			sig_f = finufft3d2( ...
				plan.fourier_pts(1,:), ...
				plan.fourier_pts(2,:), ...
				plan.fourier_pts(3,:), ...
				-1, epsilon, ...
				sig);
		end
	end

	sig_f = cast(sig_f, precision);
end
