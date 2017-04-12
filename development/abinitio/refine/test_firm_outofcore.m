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
q=qrand(Nprojs);  % Generate random orientations for the images
log_message('Generating %d projections of size %dx%d',Nprojs,n,n);
projs=cryo_project(map,q,n);  % Generate projections of map
projs=permute(projs,[2,1,3]); % For backward compatability. Will be removed 
                              % in the future.
log_message('Projections generated successully');

%% Compute rotations of the simulated projections
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end


%% Reconstrut in-core
log_message('Reconstructing density in-core...');
[ v1, v_b1, kernel1,err1,iter1,flag1] = recon3d_firm( projs,trueRs,[], 1e-6, 30, zeros(n,n,n));
log_message('Finished reconstructing density.');
ii=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative magnitude of imaginary components = % e',ii);
v1=real(v1);

%% Reconstruct out-of-core
log_message('Reconstructing density out-of-core...');
projs_fname='fname.mrc';
mat2mrc(projs_fname,projs);
[ v2, v_b2, kernel2,err2,iter2,flag2] = recon3d_firm_outofcore( projs_fname,trueRs,[], 1e-6, 30, zeros(n,n,n));
log_message('Finished reconstructing density.');
ii=norm(imag(v2(:)))/norm(v2(:));
log_message('Relative magnitude of imaginary components = % e',ii);
v2=real(v2);

%% Compare results
err=norm(v1(:)-v2(:))/norm(v1(:));
fprintf('Difference in reconstructed volume = %e \t',err);
if err<5*eps('single')
    fprintf('OK\n');
else
    fprintd('FAILED');
end

err=norm(v_b1(:)-v_b2(:))/norm(v_b1(:));
fprintf('Difference in kernels = %e \t',err);
if err<5*eps('single')
    fprintf('OK\n');
else
    fprintd('FAILED');
end

err=norm(err1(:)-err2(:))/norm(err1(:));
fprintf('Difference in residuals = %e \t',err);
if err<5*eps('single')
    fprintf('OK\n');
else
    fprintd('FAILED');
end
