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

%% Initialize parameters for simulation
% The next 3 parameters are used to construct a discretezation of SO(3)
% which will be used to estimate the relative roatitions between each pair 
% of images. 
grid_res = 1200;  % Number of points to sample on the sphere as projection 
                  % directions. Default = 1200. 
eq_min_dist = 15;  % Defines which projection directions are considered to
                  % be 'Equator directions'and are to be filtered. Default = 15. 
inplane_res = 5;  % sampling resolution of angles of in plane rotations. 
                  % Default = 5. (i.e. if ntheta == 360, then sample 72 
                  % angulary equispaced inplane rotations for each 
                  % projection direction). 
ntheta = 360;     % Angular resolution (number of Fourier rays) used when 
                  % searching for common lines between pairs of images.
%% Run algorithm
results = runD2(projs,gpuIdx,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,ntheta,doFilter,s);

