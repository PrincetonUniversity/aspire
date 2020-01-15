
initpath; % Need to input local directory of aspire inside initpath file!!! 
delete(gcp('nocreate'));

%% Volume for simulation
res = 89;
vol = load('/home/eitanr/cryo/D2ForAspire/7770.mat');
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
nproj = 100;
[projs,Rijs_gt,q,ref_shifts]=genDataForSimulation(vol,...
    nproj,max_shift,1,snr,s,0);

doFilter=1;
gpuIdx = 1:2;
nCpu = maxNumCompThreads;

pixA=1.896;
cutoff=0.143;
debugParam=struct('q',q,'vol',vol,'pixA',1.896,'cutoff',0.143);

grid_res = 1200;
eq_min_dist = 7;
inplane_res = 5;
ntheta = 360;

[results]=runD2(projs,gpuIdx,nCpu,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,ntheta,doFilter,s,Rijs_gt,q,debugParam);
