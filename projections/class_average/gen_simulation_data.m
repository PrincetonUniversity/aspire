%gen_simulation_data
%Generate 10^4 clean projection images from 70S ribosome volume (in ./simulation/volume.mat). Image size
%is 129x129 pixels.

fname = mfilename('fullpath');
[pathstr,~,~] = fileparts(fname); 

pathstr=fullfile(pathstr,'simulation');
if ~exist(pathstr,'dir')
    if ~mkdir(pathstr)
        error('Cannot create directory %s',pathstr);
    end
end

K = 10000; %K is the number of images

% Set parameters
tol=1.0e-5; % Required accuracy
batch_size=1000;
load volume.mat

% Estimate timing for generating images
log_message('Estimating time required for generating images');
K1=100;
rots = rand_rots(K1);
tic; cryo_project(vol,rots,size(vol,1),tol,10); t1=toc;
K2=200;
rots = rand_rots(K2);
tic; cryo_project(vol,rots,size(vol,1),tol,10); t2=toc;

% Extrapolate to K images
est_time=(K-K2)/(K1-K2)*t1+(K-K1)/(K2-K1)*t2;
log_message('Estimated for for generating %d images is %3.0f seconds',K,est_time);
log_message('Generating images...');

max_shift = 0; %maximum shift range
step_size = 1;
initstate;
rots = rand_rots(K);
shifts=round((rand(K,2)-1/2)*2*max_shift/step_size)*step_size;
tic
projections=cryo_project(vol,rots,size(vol,1),1.0e-5,1000);
t=toc;
log_message('Finished in %1.0f seconds',t);
save(fullfile(pathstr,'clean_data'), '-v7.3', 'projections', 'rots', 'shifts')
clear

