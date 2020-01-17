%% Initialize path
delete(gcp('nocreate'));

%% Initialize images and ML data
mapname='/scratch/yoelsh/tmp/stack.mrcs';
projstack=ReadMRC(mapname);
sidx=1:size(projstack,3);%9002:2:10000;
projs=projstack(:,:,sidx);
nr=size(projs,1);
max_shift_ratio=0.15;
max_shift=round(nr*max_shift_ratio);
shift_step=1;
vol=ReadMRC(mapname); %Input vol from mrc for comparison
s=rng();

%% Generate lookup Data
grid_res=1200;
eq_min_dist=15;
inplane_res=5;

%% Initialize parameters for algorithm and run
ntheta=360;
doFilter=1;
gpuIdx = 1:2;
nCpu = maxNumCompThreads;

%Run algorithm
[results]=runD2(projs,gpuIdx,nCpu,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,ntheta,doFilter,s);

