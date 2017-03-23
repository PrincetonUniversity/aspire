% NUFFT_SET_POINTS Set Fourier points of an NUFFT plan
%
% Usage
%    plan = nufft_set_points(plan, fourier_pts);
%
% Input
%    plan: An NUFFT plan.
%    fourier_pts: The points in Fourier space where the Fourier transform is to
%       be calculated, arranged as an N-by-d array, where d is the dimension-
%       ality of the plan. These need to be in the range [-pi, pi] in each dimen-
%       sion.

function plan = nufft_set_points(plan, fourier_pts)
	dims = numel(plan.sz);

	if ~all(size(fourier_pts)==[plan.num_pts dims])
		error('Frequencies ''fourier_pts'' array must be of the size N-by-d.');
	end

	plan.fourier_pts = fourier_pts;

	if plan.lib_code == 3
		nfft_set_x(plan.nfft_plan_id, 1/(2*pi)*fourier_pts');
		nfft_precompute_psi(plan.nfft_plan_id);
	end
end
