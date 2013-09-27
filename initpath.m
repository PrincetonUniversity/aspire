% Setup path

fname = mfilename('fullpath');
[pathstr,~,~] = fileparts(fname); % Find where the package installed.

addpath(genpath(fullfile(pathstr,'abinitio')));
addpath(genpath(fullfile(pathstr,'common')));
addpath(genpath(fullfile(pathstr,'examples')));
addpath(genpath(fullfile(pathstr,'fourier')));
addpath(genpath(fullfile(pathstr,'projections')))
addpath(genpath(fullfile(pathstr,'sinograms')))
addpath(genpath(fullfile(pathstr,'reconstruction')))


addpath(fullfile(pathstr,'extern','SDPLR-1.03-beta'))
run(fullfile(pathstr,'extern','irt','setup.m'))
run(fullfile(pathstr,'extern','cvx','cvx_startup.m'));
