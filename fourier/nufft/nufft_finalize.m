% NUFFT_FINALIZE Finalize NUFFT plan
%
% Usage
%    nufft_finalize(plan);
%
% Input
%    plan: An NUFFT plan.

function nufft_finalize(plan)
	if plan.lib_code == 3
		nfft_finalize(plan.nfft_plan_id);
	end
end
