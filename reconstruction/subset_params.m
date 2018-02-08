% SUBSET_PARAMS Extract subset of parameters
%
% Usage
%    params = subset_params(params, ind);
%
% Input
%    params: The imaging parameters structure from which we want to extract a
%       subset.
%    ind: A logical array or a set of indices defining the subset.
%
% Output
%    params: The parameters corresponding to the desired subset.

function params = subset_params(params, ind)
    if nargin < 2
        error('Input `ind` must be specified.');
    end

    check_imaging_params(params);

    if ndims(ind) ~= 2 || size(ind, 1) ~= 1
        error('Input `ind` must be of size 1-by-m.');
    end

    if islogical(ind)
        ind = find(ind);
    end

    if any(floor(ind) ~= ind)
        error('Input `ind` must be a logical array or an array of indices.');
    end

    params.rot_matrices = params.rot_matrices(:,:,ind);
    params.ctf_idx = params.ctf_idx(:,ind);
    params.ampl = params.ampl(:,ind);
    params.shifts = params.shifts(:,ind);
end
