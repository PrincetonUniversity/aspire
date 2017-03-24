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
	dims = numel(plan.sz);

	if ~all(size(sig_f)==[plan.num_pts 1])
		error('Input ''sig_f'' must have plan.num_pts elements.');
	end

	precision = class(sig_f);

	if plan.lib_code == 1
		if dims == 1
			sig = anudft1(sig_f, plan.fourier_pts, plan.sz);
		elseif dims == 2
			sig = anudft2(sig_f, plan.fourier_pts, plan.sz);
		elseif dims == 3
			sig = anudft3(sig_f, plan.fourier_pts, plan.sz);
		end
	elseif plan.lib_code == 2
		epsilon = 1e-10;

		sig_f = double(sig_f(:));

		if dims == 1
			sig = plan.num_pts*nufft1d1(plan.num_pts, ...
				plan.fourier_pts(1,:), ...
				sig_f, ...
				1, epsilon, plan.sz(1));
		elseif dims == 2
			sig = plan.num_pts*nufft2d1(plan.num_pts, ...
				plan.fourier_pts(1,:), ...
				plan.fourier_pts(2,:), ...
				sig_f, ...
				1, epsilon, ...
				plan.sz(1), plan.sz(2));
		elseif dims == 3
			sig = plan.num_pts*nufft3d1(plan.num_pts, ...
				plan.fourier_pts(1,:), ...
				plan.fourier_pts(2,:), ...
				plan.fourier_pts(3,:), ...
				sig_f, ...
				1, epsilon, ...
				plan.sz(1), plan.sz(2), plan.sz(3));
		end

		sig = reshape(sig, [plan.sz 1]);
	elseif plan.lib_code == 3
		sig_f = double(sig_f);
		sig_f = sig_f(:);

		nfft_set_f(plan.nfft_plan_id, sig_f);
		nfft_adjoint(plan.nfft_plan_id);
		sig = nfft_get_f_hat(plan.nfft_plan_id);

		if dims == 2
			sig = permute(reshape(sig, plan.sz), [2 1]);
		elseif dims == 3
			sig = permute(reshape(sig, plan.sz), [3 2 1]);
		end
	end

	sig = cast(sig, precision);
end
