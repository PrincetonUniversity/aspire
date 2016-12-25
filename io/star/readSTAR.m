function datablocks=readSTAR(fname)
% READSTAR  Read data from a file in STAR format.
% 
% datablocks=readSTAR(fname)
% datablocks=readSTAR(fname,nrecords)
%
% Read the file fname, given in STAR format, into the structure datablocks.
% The returned variable datablocks is an array, with a structure for each
% data block in the file. If there is only one data block, then datablocks
% is a structure. Each data block structure has the following fields:
%   1. name     The name of the data block. A string of the form
%               'data_XXX'.
%   2. labels   A cell array describing the columns in the data table. In
%               the data table, each row corresponds to a data record. The
%               number of columns in each row is equal to the number of
%               labels. The comments in the labels (designated by '#') are
%               kept.
%   3. data     A cell array whose numer of elements is equal to the number
%               of data record. data{k} is a struct whose fields are the
%               labels (without the comments).
%
% Example:
%   datablocks=readSTAR('particles.star');
%
% Yoel Shkolnisky, July 2014.

blockidx=0; % Index of the current data block
datablocks=struct([]);
loopstate=0; % Are we currently reading a data table?
             % loopstate: 0 Not in loop
             %            1 Reading labels
             %            2 Reading data
lineno=1;    % Current line number in the STAR file.
dataidx=1;   % How many data rows were read.

% Get size of fname
finfo = dir(fname);

if isempty(finfo)
    error('Cannot find %s',fname);
end
fsize = finfo.bytes;
bytesread=0;

fid=fopen(fname,'r');

if fid<0
    error('Failed to open %s',fname);
end

printProgressBarHeader;

tline = fgets(fid); % Read a line form the file
while ischar(tline)
    bytesread=bytesread+numel(tline);
%     if mod(dataidx,1000)==0
%         disp(dataidx);
%     end
    progressTicFor(bytesread,fsize);
    
    tline=strtrim(tline);
    if ~isempty(tline) % if empty line. Skip to the next one
        if strfind(tline,'data_')
            % Beginning of a data block.
            % Each file must have one or more data blocks. The start of a
            % data block is defined by the keyword "data_" followed by an
            % optional string for identification (e.g., "data_images").
            blockidx=blockidx+1;
            datablocks(blockidx).name=tline;            
        elseif strfind(tline,'loop_')
            % Multiple values associated with one or more labels in a data
            % block are arranged in a table using the keyword "loop_"
            % followed by the list of labels and columns of values. The
            % values are delimited by whitespace (i.e., blanks, tabs,
            % end-of-lines and carriage returns). The loop must be followed
            % by an empty line to indicate its end. 
            loopstate=1;
            labelsidx=1; % Used to read the columns of the current data block.
            dataidx=1;   % How many data rows were read.
            datablocks(blockidx).labels=cell(1);
            datablocks(blockidx).data=cell(1);
        elseif (loopstate==1) && any(strfind(tline,'_')==1)
            % Read labels.
            % Label names always start with an underscore ("_"). Each
            % label may only be used once within each data block.
            datablocks(blockidx).labels{labelsidx}=tline(2:end); % remove underscore.
            labelsidx=labelsidx+1;
        elseif (loopstate==1) && ~any(strfind(tline,'_')==1)
            % Done reading labels, switch to reading data.
            % Data items or values can be numeric or strings of characters.
            % A string is interpreted as a single item when it doesn't
            % contain spaces.
            loopstate=2;
            labelsidx=labelsidx-1;
            
            rowdata=strsplit(tline); % Split the line at white spaces.
            if numel(rowdata)~=labelsidx
                % Check that the line has number of columns equal to the
                % number of labels.
                error('Malformed data line on line %d',lineno);
            end
            
            s=row2struct(rowdata,datablocks.labels);
            datablocks(blockidx).data{1}=s;            
            dataidx=2;
        elseif (loopstate==2)
            rowdata=strsplit(tline);            
            if numel(rowdata)~=labelsidx
                error('Malformed data line on line %d',lineno);
            end

            s=row2struct(rowdata,datablocks.labels);
            datablocks(blockidx).data{dataidx}=s;
            dataidx=dataidx+1;
        end
    else
        if (loopstate==2)
            % End of data block.
            loopstate=0;             
        end
    end
    tline = fgets(fid);
    lineno=lineno+1;
end

fclose(fid);

function s=row2struct(rowdata,labels)
% Convert a data row into a sturct, with field names given by labels.
for k=1:numel(labels)
    %[fieldname,~]=strtok(labels{k},'#');
    % strtok is very slow. So switch to strfind.
    ii=strfind(labels{k},'#');
    
    if isempty(ii) % No comment marker at end of label, so no need to strip it.
        ii=numel(labels{k})+1;
    end
    fieldname=labels{k}(1:ii-1);
    fieldname=strtrim(fieldname);
    % str2double is slow. Use sscanf
    %val=str2double(rowdata{k});
    [val,~,errmsg]=sscanf(rowdata{k}, '%g');
    if isempty(errmsg) 
        s.(fieldname) = val;
    else
        s.(fieldname) = rowdata{k};
    end
end