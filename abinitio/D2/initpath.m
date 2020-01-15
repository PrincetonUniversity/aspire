% Setup path

fname =  '/home/eitanr/cryo/aspire/'; %;mfilename('fullpath');
[pathstr,~,~] = fileparts(fname); % Find where the package installed.

addpath(genpath(fullfile(pathstr,'abinitio')));
addpath(genpath(fullfile(pathstr,'common')));
addpath(genpath(fullfile(pathstr,'examples')));
addpath(genpath(fullfile(pathstr,'fourier')));
addpath(genpath(fullfile(pathstr,'io')));
addpath(genpath(fullfile(pathstr,'projections')));
addpath(genpath(fullfile(pathstr,'sinograms')));
addpath(genpath(fullfile(pathstr,'reconstruction')))
addpath(genpath(fullfile(pathstr,'refinement')));
addpath(genpath(fullfile(pathstr,'workflow')));
%addpath(genpath(fullfile(pathstr,'extern/nufftall-1.33'))); %Also iside Aspire local path. 
addpath(genpath(fullfile(pathstr,'basis')));
addpath(genpath(fullfile(pathstr,'projections/simulation')));

fname = '/home/yoel/data/work/aspire/';
[pathstr,~,~] = fileparts(fname);
addpath(genpath(fullfile(pathstr,'extern/nufftall-1.33'))); %Also iside Aspire local path. 
%addpath(fullfile(pathstr,'extern','SDPLR-1.03-beta'));
%run(fullfile(pathstr,'extern','irt','setup.m'))
clear fname pathstr