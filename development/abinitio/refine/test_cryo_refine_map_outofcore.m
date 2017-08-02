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
Nprojs=100;  % Number of projections to generate
rots=rand_rots(Nprojs);  % Generate random orientations for the images
log_message('Generating %d projections of size %dx%d',Nprojs,n_orig,n_orig);
projs=cryo_project(map,rots,n_orig);  % Generate projections of map
projs=permute(projs,[2,1,3]); % For backward compatability. Will be removed 
                              % in the future.
[projs,ref_shifts]=cryo_addshifts(projs,[],2,1);
snr=1;
projs=cryo_addnoise(projs,snr,'gaussian');
                              
log_message('Projections generated successully');

%% Compute rotations of the simulated projections
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end


%% Prepare data files
vol_fname=tempmrcname;
WriteMRC(map_downsampled,1,vol_fname);
projs_fname=tempmrcname;
WriteMRC(projs,1,projs_fname);

maxiter=5;

%% Refine in-core
map_out_step1=fullfile(tempmrcdir,'incore_step1');
map_out_step2=fullfile(tempmrcdir,'incore_step2');
mat_out=fullfile(tempmrcdir,'incore_out');
cryo_refine_map(projs_fname,vol_fname,map_out_step1,map_out_step2,mat_out,maxiter)

%% Refined out-of-core
map_out_step1=fullfile(tempmrcdir,'outofcore_step1');
map_out_step2=fullfile(tempmrcdir,'outofcore_step2');
mat_out=fullfile(tempmrcdir,'outcore_out');
cryo_refine_map_outofcore(projs_fname,vol_fname,map_out_step1,map_out_step2,mat_out,maxiter)

%% Compare result
currentsilentmode=log_silent(1);
match=1;
for k=1:maxiter
    vol_incore=ReadMRC(fullfile(tempmrcdir,sprintf('incore_step2_%d.mrc',k)));
    vol_outofcore=ReadMRC(fullfile(tempmrcdir,sprintf('outofcore_step2_%d.mrc',k)));
    err=norm(vol_incore(:)-vol_outofcore(:))/norm(vol_incore(:));
    if err>5*eps('single')
        match=0;
    end
    fprintf('iter %d err=%e\n',k,err);
end
if match
    fprintf('Test ok\n');
else
    fprintf('Test FAILED\n');
end
log_silent(currentsilentmode);