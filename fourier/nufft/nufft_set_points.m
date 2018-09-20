% NUFFT_SET_POINTS Set Fourier points of an NUFFT plan
%
% Usage
%    plan = nufft_set_points(plan, fourier_pts);
%
% Input
%    plan: An NUFFT plan.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as an N-by-d array, where d is the dimension-
%       ality of the plan. These need to be in the range [-pi, pi] in each
%       dimension.

function plan = nufft_set_points(plan, fourier_pts)
	dims = numel(plan.sz);

	if ndims(fourier_pts) ~= 2 || any(size(fourier_pts) ~= [dims plan.num_pts])
		error('Frequencies ''fourier_pts'' array must be of the size d-by-N.');
	end

	if ismember(plan.lib_code, [1 2 3])
		plan.fourier_pts = fourier_pts;
	end

	if ismember(plan.lib_code, [4 5])
		plan.fourier_pts = mod(fourier_pts+pi, 2*pi)-pi;
	end

	if plan.lib_code == 3
		if plan.num_threads ~= 0
			orig_num_threads = omp_get_num_threads();
			omp_set_num_threads(plan.num_threads);
		end

		nfft_set_x(plan.nfft_plan_id, 1/(2*pi)*fourier_pts);
		nfft_precompute_psi(plan.nfft_plan_id);

		if plan.num_threads ~= 0
			omp_set_num_threads(orig_num_threads);
		end
	end
end
