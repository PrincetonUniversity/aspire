function [coords types]=ReadPDBAtoms(filename, hetatms)
% function [coords types]=ReadPDBAtoms(filename, hetatms);
% Fast pdb coordinate reader.  Reads the coordinates of all ATOM records or
% both ATOM and HETATM records (if hetatms=1, default).  Given na atoms,
% coords is a 3 x na array, and types is an na x 4 character array, with
% each type padded with trailing blanks.
% Code is modified from pdbread() in the Bioinformatics Toolbox.
% fs 26 Apr 11

if nargin<2
    hetatms=1;
end;
if hetatms>0
    atomText={'ATOM','HETATM'};
else
    atomText='ATOM';
end;

fid = fopen(filename,'r');

if fid == -1,
    error('ReadPDBAtoms:CouldNotOpenFile',...
        'Could not open file %s.', filename);
end;

theLines = textscan(fid,'%s','delimiter','\n');
theLines = theLines{1};  % textscan returns a 1x1 cell containing a cell array.
fclose(fid);

% %%
% % Remove the lines corresponding to other models when requested a
% % particular model
% if modelNum ~= 0
%     modelStarts = find(strncmp(theLines,'MODEL ',6));
%     modelEnds = find(strncmp(theLines,'ENDMDL ',7));
%     modelMatch = ['^MODEL \s*' num2str(modelNum) '\s'];
%     h = find(~cellfun(@isempty,regexp(theLines(modelStarts),modelMatch,'once')));
%     if isempty(h)
%        warning('Bioinfo:pdbread:NoModelEntry',...
%                'PDB file does not contain Model %d. Reading the entire file.', modelNum);
%     else
%        q = [1:modelStarts(1)-1 modelStarts(h):modelEnds(h) modelEnds(end)+1:numel(theLines)];
%        theLines = theLines(q);
%     end
% end
nl=numel(theLines);

% Count the number of atom coordinates
na=0;
for i=1:nl
    tline = theLines{i};
    if ischar(tline)
        len = length(tline);
        if len>5 && any(strcmpi(deblank(tline(1:6)),atomText))
            na=na+1;
        end;
    end;
end;


coords=zeros(3,na);
types=char(zeros(na,4));
na=0;
for i=1:nl
    tline0 = theLines{i};
    len = length(tline0);
    
    %     % RCSB web site format requires each line to have 80 characters. This
    %     % avoids exceeding the matrix dimension for lines with
    %     % Less than 80 characters.
    tline = [tline0 blanks(80-len)];
    Record_name = deblank(tline(1:6));
    
    switch upper(Record_name)
        %
        case atomText
            na=na+1;
            s=strtrim(tline(13:16));  % atom type
            sl=length(s);
            types(na,:)=[s blanks(4-sl)];
            coords(1,na) = sscanf(tline(31:38),'%f');
            coords(2,na) = sscanf(tline(39:46),'%f');
            coords(3,na) = sscanf(tline(47:54),'%f');
    end;
end;

% plot3(coords(1,:),coords(2,:),coords(3,:),'k.');