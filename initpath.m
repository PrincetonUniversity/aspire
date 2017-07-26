% INITPATH Initialize paths for ASPIRE toolbox
%
% Usage
%    initpath();

function initpath()
    fname = mfilename('fullpath');
    [pathstr,~,~] = fileparts(fname); % Find where the package installed.

    addpath(pathstr);

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

    if exist(fullfile(pathstr, 'extern', 'nufftall-1.33'))
        addpath(fullfile(pathstr, 'extern', 'nufftall-1.33'));
    end

    if exist(fullfile(pathstr, 'extern', 'nfft'))
        addpath(fullfile(pathstr, 'extern', 'nfft', 'lib'));
        addpath(fullfile(pathstr, 'extern', 'nfft', 'share', 'nfft', ...
            'matlab', 'nfft'));
    end

    if exist(fullfile(pathstr, 'extern', 'SDPLR-1.03-beta'))
        addpath(fullfile(pathstr,'extern','SDPLR-1.03-beta'));
    end
end
