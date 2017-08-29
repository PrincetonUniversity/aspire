function [projections, noisy_projections, shifts, rots] = ...
    cryo_gen_projections(n,K,SNR,max_shift,shift_step,ref_shifts,rots_ref,precision)
%
% Simulate projections. 
% Same functionality as gen_projections_v2, but uses the new simulation
% code.
% See test_cryo_gen_projections.m for an example.
%
% Input parameters:
%   n       Size of of each projection (nxn).
%   K       Number of projections to generate.
%   SNR     Signal-to-noise-ratio of each projection, defined as the
%           variance of the clean projection divded by the variance of the
%           noise. 
%   max_shift Maximal random shift (in pixels) introduced to each
%           projection. max_shift must be an integer and the resulting
%           random shifts are integers between -max_shift to +max_shift.
%           Default is 0 (no shift). Pass [] is ref_shifts are given.
%   shift_step Resolution used to generate shifts. step_size=1 allows
%           for all integer shifts from -max_shift to max_shift (in both x
%           and y directions). step_size=2 generates shifts between
%           -max_shift and max_shift with steps of 2 pixels, that is,
%           shifts of 0 pixels, 2 pixels, 4 pixels, and so on. Default
%           is 0. Pass [] is ref_shifts are given.
%   ref_shifts  A two column table with the 2D shift of each projection.
%           Number of rows must be equal to the number of proejctions. If
%           this parameter is provided, max_shift and shift_step are
%           ignored.
%   rots_ref  An array of rotation matrices used to generate the projections.
%           The size of the third dimension must equal the number of
%           projections. If this parameter is missing, the functions draws
%           rotation matrices uniformly at random.
%   precision   Accuracy of the projections. 'single' or 'double'. Default
%           is 'single' (faster)
%
% Output parameters:
%   projections         Clean (shifted) simulated projections.
%   noisy_projections   Noisy (shifted) simulated projections.
%   shifts              A two column table with the 2D shift introduced to
%            each projections. If ref_shift is given, this is equal to
%            ref_shifts, otherwise, the shifts are random.
%   rots     A 3-by-3-by-K array of rotation matrices used to generate each of
%            the projections. If rots_ref is given, this is equal to that array.
%            Otherwise, the rotations re random and uniformly distributed on
%            SO(3).
%
% Yoel Shkolnisky, September 2013.

if ~exist('precision','var')
    precision='single';
end

if exist('rots_ref','var')
    if ~isempty(rots_ref)
        rots = rots_ref;
    end
else
    %initstate; % So we get the same results every time for reproducibility.
    log_message('For reproducible results, call initstate before calling cryo_gen_projections ');
    rots = rand_rots(K);
end

shifts=[];
if exist('ref_shifts','var')
    if ~isempty(ref_shifts)
        shifts=ref_shifts;
    end
end

if ~exist('max_shift','var')
    max_shift=0;
end

if ~exist('shift_step','var')
    shift_step=0;
end

% Generate clean projections
load cleanrib
log_message('Generating clean projections');
projections=cryo_project(volref,rots,n,precision);


projections=permute(projections,[2 1 3]); % Swap dimensions for compitability with old gen_projections.

% Add shifts
if ~isempty(shifts)
    log_message('Adding user-provided shifts to projections');
    projections=cryo_addshifts(projections,shifts);
else
    log_message('Adding randomly-generated shifts to projections');
    [projections,shifts]=cryo_addshifts(projections,[],max_shift,shift_step);
end

% Add noise
log_message('Adding noise to projections');
noisy_projections=cryo_addnoise(projections,SNR,'gaussian');
