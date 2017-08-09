% INITPATH Initialize paths for ASPIRE toolbox
%
% Usage
%    initpath();

% Find where the package installed.
[__pathstr, ~, ~] = fileparts(mfilename('fullpath'));

addpath(__pathstr);

addpath(genpath(fullfile(__pathstr,'abinitio')));
addpath(genpath(fullfile(__pathstr,'basis')));
addpath(genpath(fullfile(__pathstr,'common')));
addpath(genpath(fullfile(__pathstr,'examples')));
addpath(genpath(fullfile(__pathstr,'fourier')));
addpath(genpath(fullfile(__pathstr,'io')))
addpath(genpath(fullfile(__pathstr,'projections')))
addpath(genpath(fullfile(__pathstr,'sinograms')))
addpath(genpath(fullfile(__pathstr,'reconstruction')))
addpath(genpath(fullfile(__pathstr,'refinement')))
addpath(genpath(fullfile(__pathstr,'workflow')))

if exist(fullfile(__pathstr, 'extern', 'nufftall-1.33'))
    addpath(fullfile(__pathstr, 'extern', 'nufftall-1.33'));
end

if exist(fullfile(__pathstr, 'extern', 'nfft'))
    addpath(fullfile(__pathstr, 'extern', 'nfft', 'lib'));
    addpath(fullfile(__pathstr, 'extern', 'nfft', 'share', 'nfft', ...
        'matlab', 'nfft'));
end

if exist(fullfile(__pathstr, 'extern', 'SDPLR-1.03-beta'))
    addpath(fullfile(__pathstr,'extern','SDPLR-1.03-beta'));
end

clear __pathstr;
