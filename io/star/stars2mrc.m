function stars2mrc(starwd,starname)
%
% STARS2MRC     Combine multiple STAR file into a single file.
%
% stars2mrc(starwd,starname)
%   Create an STAR file named starname by combining multiple STAR files,
%   whose names are specified as a single string (using wildcards) 
%   by starwd. 
%
% Example:
%   stars2mrc('shiny_*.star','combined.star');
%
% Yoel Shkolnisky, July 2017.

if ~exist('starname','var')
    error('Name of combined STAR not given');
end


% First pass, create combined STAR file
starfiles=dir(starwd);

for fidx=1:numel(starfiles)
    currstarfile=fullfile(starfiles(fidx).folder,starfiles(fidx).name);
    log_message('Loading %s',currstarfile);
    datablock=readSTAR(currstarfile);
    
    if fidx==1
        % First STAR file. Use its header to create the output STAR file.
        combineblocks=createSTARdata(-1,datablock.labels{:});
    end
    N=numel(datablock.data); % Number of images to read

    % Write STAR data into combined file
    log_message('Copying STAR data into %s',starname);
    printProgressBarHeader;
    for i=1:N
        progressTicFor(i,N);        
        % Update STAR records
        kv=struct2keyvalue(datablock.data{i});
        combineblocks=addrecordtoSTARdata(combineblocks,-1,kv{:});
    end
end
% Saving combined star data
log_message('Saving %s',starname);
writeSTAR(combineblocks,starname);
