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
%run(fullfile(pathstr,'extern','irt','setup.m'))
%run(fullfile(pathstr,'extern','cvx','cvx_startup.m'));

if exist(fullfile(pathstr, 'extern', 'nufftall-1.33'))
    addpath(fullfile(pathstr, 'extern', 'nufftall-1.33'));
end

if exist(fullfile(pathstr, 'extern', 'nfft'))
    addpath(fullfile(pathstr, 'extern', 'nfft', 'lib'));
    addpath(fullfile(pathstr, 'extern', 'nfft', 'share', 'nfft', ...
        'matlab', 'nfft'));
end
