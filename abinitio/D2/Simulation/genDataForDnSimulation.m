
function [noisy_projs,Rijs_gt,q,ref_shifts]=genDataForSimulation(vol,...
    n_projs,max_shift,shift_step,snr,s,eqFilterAngle,shifts,grid_in,n)
%% Analyse input
if exist('eqFilterAngle','var')
    if isempty(eqFilterAngle)
        eqFilterAngle=0;
    end
else
    eqFilterAngle=0;
end
if ~exist('s','var')
    s=rng('shuffle');
end
if ~exist('snr','var')
    snr=1000000;
end
if ~exist('shift_step','var')
    shift_step=1;
end
if ~exist('max_shift','var')
    max_shift=0;
end

if isempty(gcp('nocreate'))
    parpool('local',12);
end

%% Generate projections from volume
if ~exist('grid_in','var')
    [test_data]=genRotations_Random(n_projs,eqFilterAngle,s);
else
    [eq_idx,~]=markEquatorsDn2(squeeze(grid_in(:,3,:)),eqFilterAngle,n);
    grid_in=grid_in(:,:,~eq_idx);
    nrot=size(grid_in,3);
    q=zeros(4,nrot);
    for i=1:nrot
        q(:,i)=rot_to_q(grid_in(:,:,i));
    end
    test_data=struct('rots_grid',grid_in,'grid_stats',[],'nrot',size(grid_in,3),...
    'eq_filter_angle',eqFilterAngle,'pol_rep',[],'q',q);
end

test_grid=test_data.rots_grid;
[Rijs_gt,~]=genRijsFromGridD2(test_grid);
test_grid=permute(test_grid,[2,1,3]);
%[Rijs_gt,~]=genRijsFromGridD2(permute(test_grid,[2,1,3]));
projs = cryo_project(vol,test_grid);
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
q=test_data.q;
