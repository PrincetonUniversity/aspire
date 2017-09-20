clear;

%% Read a density map from EMDB
EMDID=2275;  % Code in EMDB of the map to download
log_message('Downloading map with EMDID=%d from EMDB',EMDID);
mapname=cryo_fetch_emdID(EMDID);  % Returns full path to downloaded map
map=ReadMRC(mapname); % Read downloaded map
map=single(map);  % Convert to single precision - saves memory
log_message('Map downloaded successully');


%% Downsample map
n_orig=size(map,1);  % Size of the original map
n=65;   % Size of the resampled (downsampled map)
log_message('Resampling map from %dx%dx%d to %dx%dx%d',...
    n_orig,n_orig,n_orig,n,n,n);
map=cryo_downsample(map,[n n n]); % Downsamples map
log_message('Map resampled successully');


%% Project map
Nprojs=500;  % Number of projections to generate
rots = rand_rots(Nprojs);  % Generate random orientations for the images
log_message('Generating %d projections of size %dx%d',Nprojs,n,n);
projs=cryo_project(map,rots,n);  % Generate projections of map
projs=permute(projs,[2,1,3]); % For backward compatability. Will be removed 
                              % in the future.
log_message('Projections generated successully');

%% Compute rotations of the simulated projections
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end


%% Reconstrut
log_message('Reconstructing density...');
[ v, v_b, kernel ,err, iter, flag] = recon3d_firm( projs,trueRs,[], 1e-6, 30, zeros(n,n,n));
log_message('Finished reconstructing density.');
ii=norm(imag(v(:)))/norm(v(:));
log_message('Relative magnitude of imaginary components = % e',ii);
v=real(v);

%% Print results
log_message('Correlation of original and reconstructed volumes = %7.2f',corr(map(:),v(:)));
pixA=1.77*n_orig/n; % From emd-2275.xml
log_message('Pixel size = %d',pixA);
plotFSC(map,v,0.143,pixA);