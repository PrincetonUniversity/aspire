% Simulate projections images
%
% The script demonstrates how to generate projections images from a given density map.
%
% Follow comments below.
%
% Yoel Shkolnisky, March 2017.

%% Read a density map from EMDB
EMDID=2275;  % Code in EMDB of the map to download
log_message('Downloading map with EMDID=%d from EMDB',EMDID);
mapname=cryo_fetch_emdID(EMDID);  % Returns full path to downloaded map
map=ReadMRC(mapname); % Read downloaded map
map=single(map);  % Convert to single precision - saves memory
log_message('Map downloaded successully');


%% Downsample map
n_orig=size(map,1);  % Size of the original map
n=89;   % Size of the resampled (downsampled map)
log_message('Resampling map from %dx%dx%d to %dx%dx%d',...
    n_orig,n_orig,n_orig,n,n,n);
map=cryo_downsample(map,[n n n]); % Downsamples map
log_message('Map resampled successully');

%% Project map
Nprojs=1000;  % Number of projections to generate
rots = rand_rots(Nprojs);  % Generate random orientations for the images
log_message('Generating %d projections of size %dx%d',Nprojs,n,n);
projs=cryo_project(map,rots,n);  % Generate projections of map
projs=permute(projs,[2,1,3]); % For backward compatability. Will be removed 
                              % in the future.
log_message('Projections generated successully');                              

%% Add shifts to the projections
 log_message('Adding shifts to projections');
 max_shift=5;   % Maximal shift in each axis
 shift_step=1; 
 [projs,shifts]=cryo_addshifts(projs,[],max_shift,shift_step);
 log_message('Shifts added successully')
 log_message('Shifts added to the first 5 projections:')
 for k=1:5
     log_message('\t [%2d,%2d]',shifts(k,1),shifts(k,2));
 end
 
%% Add noise to the projections
log_message('Adding noise to projections');
SNR=1/16;
noisy_projs=cryo_addnoise(projs,SNR,'gaussian');
log_message('Noise added successully');

%% Save and display results
mrcdir=tempmrcdir;
fname_clean=fullfile(mrcdir,'clean.mrc');
fname_noisy=fullfile(mrcdir,'noisy.mrc');
log_message('Saving clean projections to %s',fname_clean);
WriteMRC(projs,1,fname_clean);
log_message('Saving noisy (and shifted) projections to %s',fname_noisy);
WriteMRC(noisy_projs,1,fname_noisy);
log_message('Stacks saved successully');

%% Viewing instructions
log_message('To load and view the stacks use:');
log_message('\t projs=ReadMRC(''%s'');',fname_clean);
log_message('\t figure; viewstack(projs,5,5);'); % Show the first 25 images on a 5x5 grid
log_message('\t projs=ReadMRC(''%s'');',fname_noisy);
log_message('\t figure; viewstack(projs,5,5);'); 

