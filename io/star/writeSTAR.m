function writeSTAR(datablock,fname)
% WRITESTAR     Write data to file in STAR format.
%
% writeSTAR(datablock,fname)
%   Write the data in datablocks into a STAR file fname. 
%
% See
% http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Conventions_%26_File_formats
% for information about STAR files.
%
% Yoel Shkolnisky, August 2015.

% Open file for writing.
if exist(fname,'file')==2
    [pathstr,name,ext]=fileparts(fname);
    name=sprintf('%s_sav',name);
    bakname=fullfile(pathstr,[name ext]);
    copyfile(fname,bakname);
    log_message('File %s alreadey exists. Moving it to %s',fname,bakname);    
end

fid=fopen(fname,'w');
if fid<0
    error('Failed to open %s. Aborting...',fname);
end
    
% The start of a data block is defined by the keyword "data_" followed by
% an optional string for identification (e.g., "data_images").
% If not name is given in the datablock, use "data_".

dbname='data_';
if isfield(datablock,'name')
    dbname=datablock.name;
end
fprintf(fid,'%s\n',dbname);
    
% Multiple values associated with one or more labels in a data block can be
% arranged in a table using the keyword "loop_" followed by the list of
% labels and columns of values. The values are delimited by whitespace
% (i.e., blanks, tabs, end-of-lines and carriage returns). The loop must be
% followed by an empty line to indicate its end.
fprintf(fid,'loop_\n');

% Label names always starts with an underscore ("_"). Each label may only
% be used once within each data block.
nlabels=numel(datablock.labels);
for k=1:nlabels;    
    label=datablock.labels{k};
    if label(1)~='_'
        label=strcat('_', label);
    end
    fprintf(fid,'%s\n',label);
end
   
% Write columns of values.
printProgressBarHeader;
ndata=numel(datablock.data);
for k=1:ndata
    progressTicFor(k,ndata);
    
    s=datablock.data{k};
    s_fieldnames=fieldnames(s);
    
    % Write current record
    for j=1:nlabels
        val=s.(s_fieldnames{j});
        if ischar(val)
            fprintf(fid,'%s',val);
        else % Val is a number
            fprintf(fid,'%.6f',val);
        end
        if j<nlabels
            fprintf(fid,'\t');
        end
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n'); % Close "loop_"

fclose(fid);