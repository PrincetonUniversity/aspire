% NUFFT_ADJOINT Apply adjoint transform according to NUFFT plan
%
% Usage
%    sig = nufft_adjoint(plan, sig_f);
%
% Input
%    plan: An NUFFT plan.
%    sig_f: A non-uniform transform of a signal.
%
% Output
%    sig: The adjoint transform of sig_f.

function sig = nufft_adjoint(plan, sig_f)
	if ~isfield(plan, 'lib_code') || ~isfield(plan, 'sz') || ...
		~isfield(plan, 'num_pts')
		error('Input ''plan'' is not a valid NUFFT plan.');
	end

	if ~isfield(plan, 'fourier_pts')
		error('Plan has not been initialized with Fourier points.');
	end

	dims = numel(plan.sz);

	if asize(sig_f, 1) ~= plan.num_pts
		error('Input ''sig_f'' must be of size plan.num_pts-by-L.');
	end

	[sig_f, sz_roll] = unroll_dim(sig_f, 2);

	epsilon = max(plan.epsilon, eps(class(sig_f)));

	if plan.lib_code == 1
		if dims == 1
			sig = anudft1(sig_f, plan.fourier_pts, plan.sz);
		elseif dims == 2
			sig = anudft2(sig_f, plan.fourier_pts, plan.sz);
		elseif dims == 3
			sig = anudft3(sig_f, plan.fourier_pts, plan.sz);
		end
	elseif plan.lib_code == 2
		sig_f = double(sig_f);

		% NUFFT errors if we give epsilon in single precision.
		epsilon = double(epsilon);

		L = size(sig_f, 2);

		sz_out = plan.sz;
		if numel(plan.sz) == 3
			sz_out = [sz_out(1) prod(sz_out(2:3))];
		end

		sig = zeros([sz_out L], class(sig_f));

		if dims == 1
			for ell = 1:L
				sig(:,ell) = plan.num_pts*nufft1d1(plan.num_pts, ...
					plan.fourier_pts(1,:), ...
					sig_f(:,ell), ...
					1, epsilon, plan.sz(1));
			end
		elseif dims == 2
			for ell = 1:L
				sig(:,:,ell) = plan.num_pts*nufft2d1(plan.num_pts, ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					sig_f(:,ell), ...
					1, epsilon, ...
					plan.sz(1), plan.sz(2));
			end
		elseif dims == 3
			for ell = 1:L
				sig(:,:,ell) = plan.num_pts*nufft3d1(plan.num_pts, ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					plan.fourier_pts(3,:), ...
					sig_f(:,ell), ...
					1, epsilon, ...
					plan.sz(1), plan.sz(2), plan.sz(3));
			end
		end

		sig = reshape(sig, [plan.sz L]);
	elseif plan.lib_code == 3
		if ~isfield(plan, 'nfft_plan_id')
			error('Input ''plan'' is not a valid NUFFT plan.');
		end

		sig_f = double(sig_f);

		L = size(sig_f, 2);

		sig = zeros(prod(plan.sz), L, class(sig_f));

		if plan.num_threads ~= 0
			orig_num_threads = omp_get_num_threads();
			omp_set_num_threads(plan.num_threads);
		end

		for ell = 1:L
			nfft_set_f(plan.nfft_plan_id, sig_f(:,ell));
			nfft_adjoint(plan.nfft_plan_id);
			sig(:,ell) = nfft_get_f_hat(plan.nfft_plan_id);
		end

		if plan.num_threads ~= 0
			omp_set_num_threads(orig_num_threads);
		end

		if dims == 2
			sig = permute(reshape(sig, [plan.sz L]), [2 1 3]);
		elseif dims == 3
			sig = permute(reshape(sig, [plan.sz L]), [3 2 1 4]);
		end
	elseif plan.lib_code == 4 || plan.lib_code == 5
		sig_f = double(sig_f);

		% FINUFFT errors if we give epsilon in single precision.
		epsilon = double(epsilon);

		if plan.num_threads ~= 0
			orig_num_threads = omp_get_num_threads();
			omp_set_num_threads(plan.num_threads);
		end

		L = size(sig_f, 2);

		sig = zeros([plan.sz L], class(sig_f));

		if dims == 1
			for ell = 1:L
				sig(:,ell) = finufft1d1( ...
					plan.fourier_pts(1,:), ...
					sig_f(:,ell), ...
					1, epsilon, plan.sz(1));
			end
		elseif dims == 2 && plan.lib_code == 4
			for ell = 1:L
				sig(:,:,ell) = finufft2d1( ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					sig_f(:,ell), ...
					1, epsilon, ...
					plan.sz(1), plan.sz(2));
			end
		elseif dims == 2 && plan.lib_code == 5
			sig = finufft2d1many( ...
				plan.fourier_pts(1,:), ...
				plan.fourier_pts(2,:), ...
				sig_f, ...
				1, epsilon, ...
				plan.sz(1), plan.sz(2));
		elseif dims == 3
			for ell = 1:L
				sig(:,:,:,ell) = finufft3d1( ...
					plan.fourier_pts(1,:), ...
					plan.fourier_pts(2,:), ...
					plan.fourier_pts(3,:), ...
					sig_f(:,ell), ...
					1, epsilon, ...
					plan.sz(1), plan.sz(2), plan.sz(3));
			end
		end

		if plan.num_threads ~= 0
			omp_set_num_threads(orig_num_threads);
		end
	end

	sig = cast(sig, class(sig_f));

	sig = roll_dim(sig, sz_roll);
end
