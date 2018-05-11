%
% log_message(varargin)
% 
% Write GCAR log message to a previosuly opened file.
%
% varargin is either in the form 'format','arg1',arg2',..., or just a
% string.
%
% Yoel Shkolnisky, November 2008
%

function log_message(varargin)
% varargin is either in the form 'format','arg1',arg2',..., or just a
% string

global log_run_silent
global log_prefix   
global log_last_prefix_time
global log_buffer
global log_max_lines;
global log_current_line
global log_bytes_written
global log_max_bytes

%%% HACK
%return;

if isempty(log_last_prefix_time)
    warning('log system was not initialized. Initalizing log to print to screen.');
    open_log(0);
end

% Do not allow log files that are too big . If more than log_max_bytes (say
% 100MB) have been written then do nothing.
if log_bytes_written>log_max_bytes
    return; 
end

newline=char(10);

if numel(log_run_silent)==0
    log_run_silent=0;
end

if numel(varargin)>1
    str=sprintf(varargin{1},varargin{2:end});
else
    str=sprintf(varargin{1});
end

% datestr is slow, so don't call it for each message, but only if the
% time since last message was more than a second. Otherwise use the old
% prefix string.
current_time=clock;
timediff=etime(current_time,log_last_prefix_time);
if timediff>1
    log_prefix=datestr(clock);
    log_last_prefix_time=current_time;
end

msg='';
if ~isempty(str) % '' stand for "print newline without a timestamp".
    msg=[log_prefix ' ' str newline];
end

if (log_current_line>log_max_lines) || (timediff>30)
    % Flush if cache full or last message was written more than 30 seconds
    % ago. Change in the future to flush log at least every 30 seconds.
    log_flush;
end

log_buffer{log_current_line}=msg;
log_current_line=log_current_line+1;
log_bytes_written=log_bytes_written+numel(msg);

if ~log_run_silent
    fprintf('%s',msg);
end
