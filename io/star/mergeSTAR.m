function mergeSTAR(outstar,inpattern,verbose)
% MERGESTAR  Combine STAR files
% 
% mergeSTAR(outstar,'instar*.mrc')
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

if ~exist('verbose','var')
    verbose=0;
end

files=dir(inpattern);
nfiles=numel(files);

if nfiles==0
    error('No files matching input pattern');
end

% Read the first file to create the output data structure
fname=fullfile(files(1).folder,files(1).name);
log_message('Reading %s to read field names',fname);
datablock=readSTAR(fname,0);

% Create a data structure with the field names into which we will insert
% the data of each record of each of the input STAR files
recfields=fieldnames(datablock.data{1});
newdatarec=cell(numel(recfields)*2,1);
for k=1:numel(recfields)
    newdatarec{2*k-1}=recfields(k);
end

nrecs=0;
% First pass - count number of records (for speedup of later updates)
log_message('First pass - counting records');
for idx=1:nfiles
    fname=fullfile(files(idx).folder,files(idx).name);
    log_message('Reading file %d/%d (%s)',idx,nfiles,fname);
    datablock=readSTAR(fname,0);
    nrecs=nrecs+numel(datablock.data);
end    

outdata=createSTARdata(nrecs,datablock.labels{:});


% Read data from all files
log_message('Second pass - reading data');
currentrec=1;
for idx=1:nfiles
    fname=fullfile(files(idx).folder,files(idx).name);
    log_message('Reading file %d/%d (%s)',idx,nfiles,fname);
    datablock=readSTAR(fname,0);
    
    % Populate values into newdatarec from all available files
    for recidx=1:numel(datablock.data)
        currentdatarec=datablock.data{recidx};
        for k=1:numel(recfields)
            newdatarec{2*k}=currentdatarec.(recfields{k});
        end
        outdata=addrecordtoSTARdata(outdata,currentrec,newdatarec{:});
        currentrec=currentrec+1;
    end

end    
log_message('Saving %s',outstar);
writeSTAR(outdata,outstar);