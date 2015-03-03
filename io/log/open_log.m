%
% open_log(logfile)
%
% Open log file for GCAR log messages
%
% Yoel Shkolnisky, November 2008
%
function open_log(logfile)
global log_fname
global log_prefix   
global log_last_prefix_time
global log_buffer
global log_max_lines;
global log_current_line;
global log_bytes_written
global log_max_bytes

if logfile==0
    log_fname=-1;  % write to screen
else    
    % Try to open the file to see everything goes fine.
    log_fname=logfile;
    logfile_fid=fopen(log_fname,'a');    
    if logfile_fid==-1
        warning('GCAR:debug','Cannot open log file');
    end
    fclose(logfile_fid);
end

log_last_prefix_time=clock;
log_prefix=datestr(log_last_prefix_time);

log_max_lines=1000;
log_current_line=1;
log_buffer=cell(log_max_lines,1);
log_bytes_written=0;
log_max_bytes=1000*2^20; % Log file should be no more than 1000MB.
