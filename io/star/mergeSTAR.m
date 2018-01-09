function mergeSTAR(outstar,inpattern)
% MERGESTAR  Combine STAR files
% 
% mergeSTAR(outstar,'in*.star')
%
% Combine multiple STAR files specified by the pattern inpattern into a
% single STAR file named outstar. Wildcards in instars are allowed.
% All input STAR files must have the same fields.
% Set verbose to nonzero to print progress messages (default 0).
%
% Example:
%   mergeSTAR('out.star','in*.star')
%
% Yoel Shkolnisky, December 2017.

files=dir(inpattern);
nfiles=numel(files);

if nfiles==0
    error('No files matching input pattern');
end

% Read the first file to create the output data structure
fname=fullfile(files(1).folder,files(1).name);
log_message('Reading %s to read field names',fname);
datablock=readSTAR(fname,0);
outdata=createSTARdata(0,datablock.labels{:});

% Read data from all files
for idx=1:nfiles
    fname=fullfile(files(idx).folder,files(idx).name);
    log_message('Reading file %d/%d (%s)',idx,nfiles,fname);
    datablock=readSTAR(fname,0);
    
    outdata.data=vertcat(outdata.data(:),datablock.data(:));
end    
log_message('Saving %s',outstar);
writeSTAR(outdata,outstar);
