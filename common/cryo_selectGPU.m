function gpuid=cryo_selectGPU(fname,verbose)
% CRYO_GETGPUCONF  Select preferred GPU (if any)
%
% gpuid=cryo_selectGPU(fname)
%   Find if any GPU is available. If not, return 0. Otherwise, if fname
%   exists, read from it the index of the preferred GPU to be used. If
%   fname not given, is empty array, or does not exist, use GPU #1.
%   If fname is a number, it is treated as a gpuid.
%   The preferred GPU (if any) is selected and initialized, so subsequent
%   functions can use the GPU.
%
% gpuid=cryo_selectGPU(fname,verbose)
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

if isempty(fname)
    gpuid=1;
    log_message('GPU detected but fname is empty. Using gpuid=1.');
elseif isscalar(fname)
    gpuid=fname;
    log_message('GPU detected and gpuid given. Using gpuid=%d.',gpuid);
elseif ~exist(fname,'file')
    log_message('GPU detected but GPU preference file not found. Using gpuid=1.');
    gpuid=1;
else
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
end

log_message('Selecting GPU#%d',gpuid);
gpuDevice(gpuid);

log_silent(currentsilentmode);

