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

	sig_sz = asize(sig);

	if ~isprefix(plan.sz, sig_sz)
		error('First dimensions of input ''sig'' must align with plan.sz.');
	end

	[sig, sz_roll] = unroll_dim(sig, dims+1);

	epsilon = max(plan.epsilon, eps(class(sig)));

	if plan.lib_code == 1
		if dims == 1
			sig_f = nudft1(sig, plan.fourier_pts);
		elseif dims == 2
			sig_f = nudft2(sig, plan.fourier_pts);
		elseif dims == 3
			sig_f = nudft3(sig, plan.fourier_pts);
		end
	elseif plan.lib_code == 2
		sig = double(reshape(sig, prod(plan.sz), []));

		% NUFFT errors if we give epsilon in single precision.
		epsilon = double(epsilon);

		L = size(sig, 2);

		sig_f = zeros(plan.num_pts, L, class(sig));

		if dims == 1
			for ell = 1:L
				sig_f(:,ell) = nufft1d2(plan.num_pts, ...
					plan.fourier_pts(1,:), ...
					-1, epsilon, plan.sz(1), sig(:,ell));
			end
		elseif dims == 2
			for ell = 1:L
				sig_f(:,ell) = nufft2d2(plan.num_pts, ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					-1, epsilon, ...
					plan.sz(1), plan.sz(2), sig(:,ell));
			end
		elseif dims == 3
			for ell = 1:L
				sig_f(:,ell) = nufft3d2(plan.num_pts, ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					plan.fourier_pts(3,:), ...
					-1, epsilon, ...
					plan.sz(1), plan.sz(2), plan.sz(3), sig(:,ell));
			end
		end
	elseif plan.lib_code == 3
		if ~isfield(plan, 'nfft_plan_id')
			error('Input ''plan'' is not a valid NUFFT plan.');
		end

		sig = double(sig);

		if dims == 2
			sig = reshape(permute(sig, [2 1 3]), prod(plan.sz), []);
		elseif dims == 3
			sig = reshape(permute(sig, [3 2 1 4]), prod(plan.sz), []);
		end

		L = size(sig, 2);

		sig_f = zeros(plan.num_pts, L, class(sig));

		for ell = 1:L
			nfft_set_f_hat(plan.nfft_plan_id, sig(:,ell));
			nfft_trafo(plan.nfft_plan_id);
			sig_f(:,ell) = nfft_get_f(plan.nfft_plan_id);
		end
	elseif plan.lib_code == 4
		sig = double(sig);

		% FINUFFT errors if we give epsilon in single precision.
		epsilon = double(epsilon);

		L = size(sig, numel(plan.sz)+1);

		sig_f = zeros(plan.num_pts, L, class(sig));

		if dims == 1
			for ell = 1:L
				sig_f(:,ell) = finufft1d2( ...
					plan.fourier_pts(1,:), ...
					-1, epsilon, sig(:,ell));
			end
		elseif dims == 2
			for ell = 1:L
				sig_f(:,ell) = finufft2d2( ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					-1, epsilon, sig(:,:,ell));
			end
		elseif dims == 3
			for ell = 1:L
				sig_f(:,ell) = finufft3d2( ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					plan.fourier_pts(3,:), ...
					-1, epsilon, ...
					sig(:,:,:,ell));
			end
		end
	end

	sig_f = cast(sig_f, class(sig));

	sig_f = roll_dim(sig_f, sz_roll);
end
