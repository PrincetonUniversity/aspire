
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

%% Generate noisy projections of a D2 symmetric volume 
max_shift_ratio=0.15; % Projection images will be shifted by up to that 
            % fraction (relative to the size of the projected volume)
max_shift=round(size(vol,1)*max_shift_ratio);
shift_step=1;
snr = 1/3; %  Signal to noise ratio of noisy simualted images. 
s = 1234; % Seed to initialize rng for reproducibility of the results.
nproj = 100; % Number of projections to generate
[projs,Rijs_gt,rots,ref_shifts]=genDataForD2Simulation(vol,...
    nproj,max_shift,1,snr,s,0);
%% Initialize parameters for simulation
doFilter=1;     % Flag to indicate whether to use a Gaussian filter to 
                % smoothen projection images. 
gpuIdx = 1:2;   % Indices of GPUs to use by the D2 algorithm.

pixA=1.896; % Pixel size of the images in the simulation. Used to compute 
            % the resolution of the reconstruced volume.
cutoff=0.143;   % Cutoff threshold for the FSC.
debugParam=struct('volref',vol,'pixA',pixA,'cutoff',cutoff);

ntheta = 360;     % Angular resolution (number of Fourier rays) used when 
                  % searching for common lines between pairs of images.

% The next 3 parameters are used to construct a discretezation of SO(3)
% which will be used to estimate the relative roatitions between each pair 
% of images. 
grid_res = 1200;  % Number of points to sample on the sphere as projection 
                  % directions. Default = 1200. 
eq_min_dist = 7;  % Defines which projection directions are considered to
                  % be 'Equator directions'and are to be filtered. Default = 7. 
inplane_res = 5;  % sampling resolution of angles of in plane rotations. 
                  % Default = 5. (i.e. if ntheta == 360, then sample 72 
                  % angulary equispaced inplane rotations for each 
                  % projection direction). 


%% Reconstruct volume from projections
results=runD2(projs,gpuIdx,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,ntheta,doFilter,s,Rijs_gt,rots,debugParam);
