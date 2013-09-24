% add paths for FIRM
fname = mfilename('fullpath');
[scriptPath,~,~] = fileparts(fname); % Find where the package installed.


addpath(scriptPath);
%fprintf('Please download the package of irt following the instruction in setup.m!\n')
% The package irt is downloaded from
% http://www.eecs.umich.edu/~fessler/code/index.html.
% This software was developed at the University of Michigan by Jeff Fessler and his students. 
% addpath 'the directory of irt'
addpath(fullfile(scriptPath,'irt'))
addpath(fullfile(scriptPath,'test'))
run(fullfile(scriptPath,'irt','setup'))
