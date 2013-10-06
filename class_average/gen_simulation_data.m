%gen_simulation_data
%Generate 10^4 clean projection images from 70S ribosome volume (in ./simulation/volume.mat). Image size
%is 129x129 pixels.

K = 10000; %K is the number of images
max_shift = 4; %maximum shift range
initstate;
q = qrand(K);
shifts = max_shift * rand(K, 2);
load volume.mat
projections=cryo_project(vol,q);
save('/simulation/clean_data', '-v7.3', 'projections', 'q', 'shifts')
clear all;


