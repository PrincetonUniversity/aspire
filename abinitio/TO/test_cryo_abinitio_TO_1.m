% Test the function cryo_abinitio_TO using simulated data with O symmetry
% 
% Yoel Shkolnisky, June 2021.

mapfile = cryo_fetch_emdID(4905); % Download 3D map file with O symmetry
volref = ReadMRC(mapfile);           % Load the map file
volref=cryo_downsample(volref,129);     

Nprojs=50;      % Due to the high order of O symmertry (24), we don't need 
                % many images. 50 images is equivalent to 50x24=1200 images.

rots=rand_rots(Nprojs); % Generate random rotations
projs=cryo_project(volref,rots);
projs=permute(projs,[2 1 3]);

projs_fname=tempmrcsname;   % Generate temporary filename to save the 
                            % projections.
vol_fname=tempmrcname;                            
params_fname=tempname;

WriteMRC(projs,512*0.65/129,projs_fname); % Save simulated projections
cryo_abinitio_TO('O',projs_fname,vol_fname,[],params_fname);

% Check reconstruction
volref=volref(1:end-1,1:end-1,1:end-1);
vol2 = ReadMRC(vol_fname);
[bestR,bestdx,vol2aligned,bestcorr]=cryo_align_vols('O',volref,vol2,1);
plotFSC(volref,vol2aligned,0.5,0.65*340/129)

% Clean up
delete(projs_fname);
delete(vol_fname);
delete(params_fname);

