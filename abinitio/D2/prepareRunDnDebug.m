
initpath; % Need to input local directory of aspire inside initpath file!!! 
delete(gcp('nocreate'));

%% Volume for simulation
res = 89;
vol = load('7770');
vol = vol.volref;
vol=cryo_downsample(vol,res);
vol=genD2fromVol(vol);
d=30;
vol2=zeros(res+d,res+d,res+d);
vol2(0.5*d+1:0.5*d+res,0.5*d+1:0.5*d+res,0.5*d+1:0.5*d+res)=vol;
vol=vol2;
clear vol2;

%% Initialize some parameters for simulation
% shift_step and max_shift are defined inside cryo_clmatrix_ML_gpu_scls.m. 
% max_shift_ratio   used to bound the the range of shifts of the images to
%                   go over when approximating common lines. 
% s         pre initialized seed for later reproducibility. 
% snr       signal to noise ratio in simualted images. 
% nproj     number of images to generate for simualtion. 
max_shift_ratio=0.15;
max_shift=round(size(vol,1)*max_shift_ratio);
shift_step=1;
snr = 1/3; 
s = rng(); % choose seed for later reproducibility. 
nproj = 200;
[projs,Rijs_gt,q,ref_shifts]=genDataForSimulation(vol,...
    nproj,max_shift,1,snr,s,0);

%% Generate lookup Data
%  genSclsScoresIdxMap_eqClass generates a lookup data for maximum
%  likelihood like scheme to estimate common lines between the images. The
%  folowing parameters are used to construct the grid: 
%  grid_res     number of sampling points on sphere for projetion directions. 
%               These are generated using the Saaf - Kuijlaars algoithm.
%               Default is 1200. 
%  inplane_res   The sampling resolution of in-plane rotations for each
%                projetion direction. Default  is 5. 
%  eq_min_dist   Width of strip around equator projection directions from
%                which we DO NOT sample directions. Default is 7. 
grid_res=1200;
eq_min_dist=7;
inplane_res=5;
lookup_data=genLookupGrid_eqClass(grid_res,eq_min_dist,inplane_res);
[scls_lookup_data]=genSelfCls(lookup_data,2);
% [oct1_ij_map,oct2_ij_map]=genSclsScoresIdxMap_eqClass(scls_lookup_data);
% scls_lookup_data.oct1_ij_map=oct1_ij_map;
% scls_lookup_data.oct2_ij_map=oct2_ij_map;
clear oct1_ij_map oct2_ij_map

%% Initialize parameters and run D2 algorithm
%   doFilter    filter images using Gaussian mask for common line
%               detection. Default is 1. 
%   gpuIdx      indices of gpu devices to use for common line detection. 
%   nCpu        size of matlab pool for parallel CPU cores. 
%   ntheta      number of fourier rays on the polar FT of each 2D
%               image. 
%   saveDir     directory to save results from algorithm. 
%   sampleName  name to be used as a prefix for results .mat files.
%   pixA        pixel size in angstrom. 
%   cutoff      Fourier shell correlations cutoff. 
%   saveIntermediate    save intermediate results from 6 stages of the
%                       algorithm. 
doFilter=1;
gpuIdx = 1:2;
nCpu = maxNumCompThreads;
ntheta=360;
saveDir=tempmrcdir;
sampleName='beta_gal_sim';
pixA=1.896;
cutoff=0.143;
saveIntermediate=1;

params=struct('max_shift_ratio',max_shift_ratio,'max_shift',max_shift,...
    'shift_step',shift_step,'doFilter',doFilter,'Rijs_gt',Rijs_gt,...,
    'gpuIdx',gpuIdx,'nCpu',nCpu,'ntheta',ntheta,...
    's',s,'q',q,'saveDir',saveDir,'sampleName',sampleName,'vol',vol,...
    'ref_shifts',ref_shifts,'pixA',pixA,'cutoff',cutoff,....
    'saveIntermediate',saveIntermediate,'scl_scores',[],'J_list_in',[]);
%Which stages to run
stages.st1=1; % Maximunm likelihood
stages.st2=1; % J-sync
stages.st3=1; % Colors sync
stages.st4=1; % Signs sync
stages.st5=1; % Reconstruction

stages.Rijs_est=[];%results.Rijs_est;
stages.Rijs_synced=[];%results.Rijs_synced;
stages.Rijs_rows=[];%results.Rijs_rows;
stages.colors=[];%results.colors;
stages.rots_est=[]; % results.rots_est;
[results]=runDn(projs,lookup_data,scls_lookup_data,params,stages);
