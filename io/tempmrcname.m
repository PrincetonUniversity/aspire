function fname=tempmrcname(isstack)
%
% TEMPMRCNAME Get teroprary filename
%
%   fname=tempmrcname
%       Returns a temporary unique filename with extension '.mrc'
%       Creates the temprary dir if needed.
%   fname=tempmrcname(isstack)
%       If isstack is nonzero, returns a filename with extension mrcs.
%
% Yoel Shkolnisky, December 2016
% Revised: Y.S. January 2018  Add isstack flag

ext='.mrc';
if ~exist('isstack','var')
    isstack=0;
end

if isstack~=0
    ext='.mrcs';
end

temporarydir=tempmrcdir;
[~,name]=fileparts(tempname);
fname=fullfile(temporarydir,[name ext]);
