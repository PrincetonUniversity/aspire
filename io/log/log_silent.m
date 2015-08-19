function silent=log_silent(gosilent)
%
% LOG_SILENT    Enable/disable log printouts to screen
%
% log_silent(gosilent)
%   Disable or disable log printouts to screen (gosilent=1 for disalble,
%   gosilent=0 for enable). Messages are still saved to log file. Returns
%   the silent mode before setting it.
%
% Yoel Shkolnisky, August 2015.

global log_run_silent

silent=log_run_silent;
log_run_silent=gosilent;