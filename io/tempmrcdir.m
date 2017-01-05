function dname=tempmrcdir
%
% TEMPMRCDIR Get teroprary directory
%
%   dname=tempmrcdir
%       Returns a temporary directory name for ASPIRE output.
%       Creates the directory if does not exist.
%
% Yoel Shkolnisky, January 2016

uname=getenv('USER');
dname=sprintf('/scratch/%s/aspire_temp',uname);
if ~exist(dname,'dir')
    log_message('Creating temporary dir %s',dname);
    mkdir(dname);
end