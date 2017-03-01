function gpuid=cryo_getGPUconf(fname,verbose)
% CRYO_GETGPUCONF  Find preferred GPU (if any)
%
% gpuid=cryo_getGPUconf(fname)
%   Find if any GPU is available. If not, return 0. Otherwise, if fname
%   exists, read from it the index of the preferred GPU to be used. If
%   fname not given, is empty array, or does not exist, use GPU #1.
%   Pass empty array for fname
%
% gpuid=cryo_getGPUconf(fname,verbose)
%    Set verbose to nonzero to print log messages.
%
% Yoel Shkolnisky, February 2017.

if ~exist('fname','var')
    fname=[];
end

if ~exist('verbose','var')
    verbose=0;
end

currentsilentmode=log_silent(verbose==0);
    
if gpuDeviceCount==0
    gpuid=0;
    log_message('No GPU found.');
    return;
end

if isempty(fname) || ~exist(fname,'file')
    gpuid=1;
    log_message('GPU detected but GPU preference file not found. Using gpuid=1.');
    return;
end

try
    fid=fopen(fname,'r');
    log_message('Opened GPU preference file %s',fname);
    gpuid=fscanf(fid,'%d');
    log_message('Using gpuid=%d',gpuid);
catch
    gpuid=1;
    log_message('Error reading gpuid from GPU preference file. Using gpuid=1');
end
fclose(fid);
    
log_silent(currentsilentmode);

