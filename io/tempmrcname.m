function fname=tempmrcname
%
% TEMPMRCNAME Get teroprary filename
%
%   fname=tempmrcname
%       Returns a temporary unique filename with extension '.mrc'
%       Creates the temprary dir if needed.
%
% Yoel Shkolnisky, December 2016

temporarydir=tempmrcdir;
[~,name]=fileparts(tempname);
fname=fullfile(temporarydir,[name '.mrc']);