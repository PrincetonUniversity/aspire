function log_flush

global log_fname
global log_buffer
global log_current_line

%%% HACK
%return;

if log_fname==-1
    return;
end
logfile_fid=fopen(log_fname,'a');
if logfile_fid==-1
    warning('GCAR:debug','Cannot open log file');
end
fprintf(logfile_fid,'%s',log_buffer{1:log_current_line-1});
fclose(logfile_fid);

log_current_line=1;




