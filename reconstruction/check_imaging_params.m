% CHECK_IMAGING_PARAMS Check imaging parameters structure
%
% Usage
%    check_imaging_params(params, L, n);
%
% Input
%    params: An imaging parameters structure containing the fields:
%          - rot_matrices: A 3-by-3-by-n array containing the rotation
%             matrices of the various projections.
%          - ctf: An L-by-L-by-K array of CTF images (centered Fourier
%             transforms of the point spread functions of the images)
%             containing the CTFs used to generate the images.
%          - ctf_idx: A vector of length n containing the indices in the
%             `ctf` array corresponding to each projection image.
%          - ampl: A vector of length n specifying the amplitude multiplier
%             of each image.
%    L: The size of the images. If left empty, this is taken from the
%       `params.ctf` field.
%    n: The number of images. If left empty, this is taken from the
%       `params.rot_matrices` field.
%
% Description
%    Checks for the necessary fields of the `params` structure and verifies the
%    internal consistency of their sizes.

function check_imaging_params(params, L, n)
    if nargin < 2
        L = [];
    end

    if nargin < 3
        n = [];
    end

    if ~isfield(params, 'rot_matrices')
        error('Parameters must have `rot_matrices` field.');
    end

    if ~isfield(params, 'ctf')
        error('Parameters must have `ctf` field.');
    end

    if ~isfield(params, 'ctf_idx')
        error('Parameters must have `ctf_idx` field.');
    end

    if ~isfield(params, 'ampl')
        error('Parameters must have `ampl` field.');
    end

    if ~isfield(params, 'shifts')
        error('Parameters must have `shifts` field.');
    end

    if isempty(L)
        L = size(params.ctf, 1);
    end

    if isempty(n)
        n = size(params.rot_matrices, 3);
    end

    if ndims(params.rot_matrices) > 3 || ...
        size(params.rot_matrices, 1) ~= 3 || size(params.rot_matrices, 2) ~= 3
        error('Field `rot_matrices` must be of size 3-by-3-by-n.');
    end

    if ndims(params.ctf) > 3 || size(params.ctf, 1) ~= L || ...
        size(params.ctf, 2) ~= L
        error('Field `ctf` must be of size L-by-L-by-K.');
    end

    if ndims(params.ctf_idx) > 2 || numel(params.ctf_idx) ~= n
       error('Field `ctf_idx` must be a vector of length n.');
    end

    if ndims(params.ampl) > 2 || numel(params.ampl) ~= n
       error('Field `ampl` must be a vector of length n.');
    end

    if any(round(params.ctf_idx) ~= params.ctf_idx) || ...
        any(params.ctf_idx < 1) || any(params.ctf_idx > size(params.ctf, 3))
        error(['Field `ctf_idx` must have positive integer entries less ' ...
            'than `size(ctf, 3)`']);
    end

    if ndims(params.shifts) > 2 || size(params.shifts, 1) ~= 2 || ...
        size(params.shifts, 2) ~= n
        error('Field `shifts` must be of the size 2-by-n.');
    end
end
