% INITPATH Initialize paths for ASPIRE toolbox
%
% Usage
%    initpath();

% Find where the package installed.
[pathstr, ~, ~] = fileparts(mfilename('fullpath'));

addpath(pathstr);

addpath(genpath(fullfile(pathstr,'abinitio')));
addpath(genpath(fullfile(pathstr,'basis')));
addpath(genpath(fullfile(pathstr,'common')));
addpath(genpath(fullfile(pathstr,'examples')));
addpath(genpath(fullfile(pathstr,'fourier')));
addpath(genpath(fullfile(pathstr,'install')));
addpath(genpath(fullfile(pathstr,'io')))
addpath(genpath(fullfile(pathstr,'projections')))
addpath(genpath(fullfile(pathstr,'sinograms')))
addpath(genpath(fullfile(pathstr,'reconstruction')))
addpath(genpath(fullfile(pathstr,'refinement')))
addpath(genpath(fullfile(pathstr,'workflow')))

if exist(fullfile(pathstr, 'extern', 'nufftall-1.33'))
    addpath(fullfile(pathstr, 'extern', 'nufftall-1.33'));
end

if exist(fullfile(pathstr, 'extern', 'nfft'))
    addpath(fullfile(pathstr, 'extern', 'nfft', 'lib'));
    addpath(fullfile(pathstr, 'extern', 'nfft', 'share', 'nfft', ...
        'matlab', 'nfft'));
end

if exist(fullfile(pathstr, 'extern', 'finufft'))
    addpath(fullfile(pathstr, 'extern', 'finufft', 'matlab'));
end

if exist(fullfile(pathstr, 'extern', 'SDPLR-1.03-beta'))
    addpath(fullfile(pathstr,'extern','SDPLR-1.03-beta'));
end

clear pathstr;
