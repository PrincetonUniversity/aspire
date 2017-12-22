% GEN_CTF_PARAMS Generate a set of CTF parameters
%
% Usage
%    ctf_params = gen_ctf_params(imaging_params, defocus);
%
% Input
%    imaging_params: A pre-filled CTF parameter structure with certain values
%       already specified, excluding the defocus field.
%    defocus: A range of defocus values in angstrom (default a uniformly
%       spaced set of 7 values from 1.5e4 to 2.5e4).
%
% Output
%    ctf_params: An array of CTF parameter structures with defocus values set
%       according to the given defocus array and the remainder set according
%       to the given imaging_params structure or filled in with the default
%       values:
%          - voltage: 200,
%          - res: 5,
%          - Cs: 2.26, and
%          - alpha: 0.07.

function ctf_params = gen_ctf_params(imaging_params, defocus)
    if nargin < 1
        imaging_params = [];
    end

    if nargin < 2 || isempty(defocus)
        defocus = linspace(1.5e4, 2.5e4, 7);
    end

    ctf_params = fill_struct(imaging_params, ...
        'voltage', 200, ...
        'resolution', 5, ...
        'defocus_u', [], ...
        'defocus_v', [], ...
        'defocus_angle', [], ...
        'Cs', 2.26, ...
        'alpha', 0.07);

    ctf_params = repmat(ctf_params, [numel(defocus) 1]);

    defocus = defocus;
    defocus_angle = zeros([numel(defocus) 1]);

    defocus = num2cell(defocus);
    defocus_angle = num2cell(defocus_angle);

    [ctf_params.defocus_u] = deal(defocus{:});
    [ctf_params.defocus_v] = deal(defocus{:});
    [ctf_params.defocus_angle] = deal(defocus_angle{:});
end
