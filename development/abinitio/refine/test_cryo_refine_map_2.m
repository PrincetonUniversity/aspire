clear;

%% Read a density map from EMDB
EMDID=2275;  % Code in EMDB of the map to download
log_message('Downloading map with EMDID=%d from EMDB',EMDID);
mapname=cryo_fetch_emdID(EMDID);  % Returns full path to downloaded map
map=ReadMRC(mapname); % Read downloaded map
map=single(map);  % Convert to single precision - saves memory
map=cryo_downsample(map,[129 129 129]);
log_message('Map downloaded successully');

if mod(size(map,1),2)==0
    map=map(1:end-1,1:end-1,1:end-1);
end

n_orig=size(map,1);  % Size of the original map

%% Downsample map
n=65;   % Size of the resampled (downsampled map)
log_message('Resampling map from %dx%dx%d to %dx%dx%d',...
    n_orig,n_orig,n_orig,n,n,n);
map_downsampled=cryo_downsample(map,[n n n]); % Downsamples map
log_message('Map resampled successully');


%% Project map
Nprojs=1000;  % Number of projections to generate
q=qrand(Nprojs);  % Generate random orientations for the images
log_message('Generating %d projections of size %dx%d',Nprojs,n_orig,n_orig);
projs=cryo_project(map,q,n_orig);  % Generate projections of map
projs=permute(projs,[2,1,3]); % For backward compatability. Will be removed 
                              % in the future.
[projs,ref_shifts]=cryo_addshifts(projs,[],2,1);
snr=1000;
projs=cryo_addnoise(projs,snr,'gaussian');
                              
log_message('Projections generated successully');

%% Compute rotations of the simulated projections
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end


%% Prepare data files
vol_fname=tempmrcname;
WriteMRC(map_downsampled,1,vol_fname);
projs_fname=tempmrcname;
WriteMRC(projs,1,projs_fname);

% Refine 
map_out_step1=fullfile(tempmrcdir,'map_out_step1');
map_out_step2=fullfile(tempmrcdir,'map_out_step2');
mat_out=fullfile(tempmrcdir,'mat_out');
cryo_refine_map(projs_fname,vol_fname,map_out_step1,map_out_step2,mat_out,30)


%% Test result
vol=ReadMRC(fullfile(tempmrcdir,'map_out_step2_1.mrc'));
plotFSC(map,vol)

