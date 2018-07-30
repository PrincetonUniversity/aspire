% INITPATH Initialize paths for ASPIRE toolbox
%
% Usage
%    initpath_development();

% Call the standard initpath of ASPIRE.
initpath;

% Add development folders
% Find where the package installed.
[pathstr, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(pathstr,'development')))
clear pathstr;
