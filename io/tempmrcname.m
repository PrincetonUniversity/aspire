function fname=tempmrcname
%
% TEMPMRCNAME Get teroprary filename
%
%   fname=tempmrcname
%       Returns a temporary unique filename with extension '.mrc'
%       Creates the temprary dir if needed.
%
% Yoel Shkolnisky, December 2016

temporarydir='/scratch/aspire_temp';
if ~exist(temporarydir,'dir')
    log_message('Creating temporary dir %s',temporarydir);
    mkdir(temporarydir);
end

[~,name]=fileparts(tempname);
fname=fullfile(temporarydir,[name '.mrc']);