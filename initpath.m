% Setup path

fname = mfilename('fullpath');
[pathstr,~,~] = fileparts(fname); % Find where the package installed.

addpath(genpath(fullfile(pathstr,'abinitio')));
addpath(genpath(fullfile(pathstr,'common')));
addpath(genpath(fullfile(pathstr,'examples')));
addpath(genpath(fullfile(pathstr,'fourier')));
addpath(genpath(fullfile(pathstr,'io')))
addpath(genpath(fullfile(pathstr,'projections')))
addpath(genpath(fullfile(pathstr,'sinograms')))
addpath(genpath(fullfile(pathstr,'reconstruction')))
addpath(genpath(fullfile(pathstr,'refinement')))
addpath(genpath(fullfile(pathstr,'workflow')))

addpath(fullfile(pathstr,'extern','SDPLR-1.03-beta'))
addpath(genpath(fullfile(pathstr,'extern','aLibs')))
run(fullfile(pathstr,'extern','irt','setup.m'))
