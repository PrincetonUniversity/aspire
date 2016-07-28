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
	dims = numel(plan.sz);

	if plan.lib_code == 1
		if dims == 1
			sig_f = nudft1(sig, plan.fourier_pts);
		elseif dims == 2
			sig_f = nudft2(sig, plan.fourier_pts);
		elseif dims == 3
			sig_f = nudft3(sig, plan.fourier_pts);
		end
	elseif plan.lib_code == 2
		epsilon = 1e-10;

		sig = double(sig(:));

		if dims == 1
			sig_f = nufft1d2(plan.num_pts, ...
				plan.fourier_pts(:,1), ...
				-1, epsilon, plan.sz(1), sig);
		elseif dims == 2
			sig_f = nufft2d2(plan.num_pts, ...
				plan.fourier_pts(:,1), ...
				plan.fourier_pts(:,2), ...
				-1, epsilon, ...
				plan.sz(1), plan.sz(2), sig);
		elseif dims == 3
			sig_f = nufft3d2(plan.num_pts, ...
				plan.fourier_pts(:,1), ...
				plan.fourier_pts(:,2), ...
				plan.fourier_pts(:,3), ...
				-1, epsilon, ...
				plan.sz(1), plan.sz(2), plan.sz(3), sig);
		end

		if isa(sig, 'single')
			sig_f = single(sig_f);
		end
	elseif plan.lib_code == 3
		sig = double(sig);

		if dims == 2
			sig = reshape(permute(sig, [2 1]), prod(plan.sz), 1);
		elseif dims == 3
			sig = reshape(permute(sig, [3 2 1]), prod(plan.sz), 1);
		end

		nfft_set_f_hat(plan.nfft_plan_id, sig);
		nfft_trafo(plan.nfft_plan_id);
		sig_f = nfft_get_f(plan.nfft_plan_id);

		if isa(sig, 'single')
			sig_f = single(sig_f);
		end
	end
end
