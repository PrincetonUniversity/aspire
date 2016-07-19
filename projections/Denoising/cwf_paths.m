fname = mfilename('fullpath');
[pathstr,~,~] = fileparts(fname); % Find where the CWF package is.
addpath(genpath(fullfile(pathstr,'cwf_functions')));
addpath(genpath(fullfile(pathstr,'cwf_scripts')));
addpath(genpath(fullfile(pathstr,'ffb')));
addpath(genpath(fullfile(pathstr,'kn_rankest')));
addpath([getenv('HOME') '/local/share/nfft/matlab/nfft']);
addpath([getenv('HOME') '/local/lib']);
