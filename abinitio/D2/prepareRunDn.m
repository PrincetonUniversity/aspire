
initpath; % Need to input local directory of aspire inside initpath file!!! 
delete(gcp('nocreate'));

%% Volume for simulation
res = 89;
D2mapfile = cryo_fetch_emdID(7770);
vol = ReadMRC(D2mapfile);
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
[projs,Rijs_gt,q,ref_shifts]=genDataForD2Simulation(vol,...
    nproj,max_shift,1,snr,s,0);

doFilter=1;
gpuIdx = 1:2;
nCpu = maxNumCompThreads;

pixA=1.896;
cutoff=0.143;
debugParam=struct('q',q,'vol',vol,'pixA',1.896,'cutoff',0.143);

[results]=runD2(projs,gpuIdx,nCpu,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,ntheta,doFilter,Rijs_gt,s,q,debugParam);
