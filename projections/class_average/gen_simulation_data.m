function gen_simulation_data
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
max_shift = 0; %maximum shift range
step_size = 1;
initstate;
rots = rand_rots(K);
shifts=round((rand(K,2)-1/2)*2*max_shift/step_size)*step_size;
load volume.mat
projections=cryo_project(vol,rots,size(vol,1),'single');
save(fullfile(pathstr,'clean_data'), '-v7.3', 'projections', 'rots', 'shifts')



