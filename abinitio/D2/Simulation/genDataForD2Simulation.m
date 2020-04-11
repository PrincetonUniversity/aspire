function [noisy_projs,Rijs_gt,rots_out,ref_shifts]=genDataForD2Simulation(vol,...
    n_projs,max_shift,shift_step,snr,s,eqFilterAngle,shifts,grid_in)
%% Generate noisy projections for D2 simulation.
%   Input Parameters: 
%   vol         .mat file with volume from which to generate 2D projection 
%               images. 
%   n_porjs     number of projections to generate. The actual number that
%               will be generated can be less if eqFilterAngle is > 0
%   max_shift   maximal spatial shift in simulated projections. Optional.
%               Default value is 0.
%   shift_step  resolution of spatial shifts. Optinal. Default value is 1.              
%   snr         signal to noise ration for simulated noise. Optional. 
%               Default value is 1000000.
%   s           seed to use for random random operations. Optional. 
%   eqFilterAngle   Angular distance to keep away from eqautor directions
%                   (3 great circles on the planes perpendicular to 
%                   symmetry axes). Optional. Default value is 1000000.
%   shifts      input shifts instead of randomly generated. Optional. 
%   grid_in     input projection directions instead of randomly generated.
%               Optinal.
%   Output: 
%   noisy_projs noisy projetion images generated from vol according to the 
%               the input parameters. 
%   rots        The rotations R1,...,Rn which were generated to simulate
%               the projection images. 
%   Rijs_gt     The realtive rotations computed from simulated rotations. 
%               Rij_gt has n 4-tuples of the form 
%               Ri^TgkRj, k=1,...,4, where gk are the group matrices of
%               D2 dihedaral symmetry (the Klein four). 
%   ref_shifts  2D array of the random shifts spatial (dx,dy) in the
%               simulated images. 


%% Default values. 
if ~exist('eqFilterAngle','var')
    eqFilterAngle=0;
    log_message('eqFilterAngle not specified. Using eqFilterAngle=%d',eqFilterAngle);
end
if ~exist('s','var')
    s=rng('shuffle');
    log_message('s not specified. Using random seed');
end
if ~exist('snr','var')
    snr=1000000;
    log_message('snr not specified. Using snr=%d',snr);
end
if ~exist('shift_step','var')
    shift_step=1;
    log_message('shift_step not specified. Using shift_step=%d',shift_step);
end
if ~exist('max_shift','var')
    max_shift=0;
    log_message('max_shift not specified. Using max_shift=%d',max_shift);
end

%% Generate projections from volume
if ~exist('grid_in','var')
    rots=genRandD2Rots(n_projs,eqFilterAngle,s);
else
    [eq_idx,~]=markEquatorsDn2(squeeze(grid_in(:,3,:)),eqFilterAngle,2);
    rots = grid_in(:,:,~eq_idx);
end

[Rijs_gt,~]=genRijsFromGridD2(rots);
rots_out = rots;
rots = permute(rots,[2,1,3]);
%[Rijs_gt,~]=genRijsFromGridD2(permute(rots,[2,1,3]));
projs = cryo_project(vol,rots);
projs = permute(projs,[2,1,3]);

%% Add noise and shifts
if ~exist('shifts','var')
    shifts=[];
end
if max_shift>0
    [shifted_projs,ref_shifts]=cryo_addshifts(projs,shifts,max_shift,shift_step);
else
    shifted_projs=projs;
    ref_shifts=[];
end
[noisy_projs,~,~,~]=...
    cryo_addnoise(shifted_projs,snr,'gaussian',s);
