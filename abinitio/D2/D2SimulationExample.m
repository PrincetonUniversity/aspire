
delete(gcp('nocreate'));

%% Generate a D2 symmetric volume for simulation
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

%% Geberate noisy projections of a D2 symmetric volume 

max_shift_ratio=0.15; % Projection images will be shifted by up to that 
            % fraction (relative to the size of the projected volume)
max_shift=round(size(vol,1)*max_shift_ratio);
shift_step=1;
snr = 1/3; %  Signal to noise ratio of noisy simualted images. 
s = 1234; % Seed for initialize rng for reproducibility of the results.
nproj = 300; % Number of projections to generate
[projs,Rijs_gt,q,ref_shifts]=genDataForD2Simulation(vol,...
    nproj,max_shift,1,snr,s,0);


%% Initialize parameters for simulation
doFilter=1;     % COMMENT
gpuIdx = 1:2;   % Indices of GPUs to use by the D2 algorithm.
nCpu = maxNumCompThreads;   % COMMENT

pixA=1.896; % Pixel size of the images in the simulation. Used to computed 
    % the resolution of the reconstruced volume.
cutoff=0.143;   % Cutoff threshold for the FSC.
debugParam=struct('q',q,'vol',vol,'pixA',pixA,'cutoff',cutoff);

grid_res = 1200;  % COMMENT
eq_min_dist = 7;  % COMMENT
inplane_res = 5;  % COMMENT
ntheta = 360;     % Angular resolution (number of Fourier ray) used when 
                  % searching for common lines between pairs of images.

%% Reconstruct volume from projections
results=runD2(projs,gpuIdx,nCpu,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,ntheta,doFilter,s,Rijs_gt,q,debugParam);
