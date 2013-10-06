function s=ReadTiffMetadata(filename)
% function s=ReadTiffMetadata(filename)
% Get the metadata from a TVIPS tiff file.  The actual image data can be
% read using imread, except that the unit16 values should be converted to
% int16.
% The returned struct has fields
% s.comment     - a string
% s.datenum     - Matlab date number
% s.datestring  - e.g. 27-Oct-2010 16:45:27
% s.volts       - e.g. 200000.
% s.pixA        - Angstroms per pixel

% Fill in default values for the returned structure
s.comment='';
s.datenum=0;
s.datestring='';
s.ecounts=0;
s.volts=0;
s.pixA=0;


f=fopen(filename,'r','ieee-le');
if f<0
    error(['The file could not be opened: ' filename])
end;

% Start reading the image file header
id=fread(f,2,'*char')';
magic=fread(f,1,'*int16');
if magic ~= 42
    error('Not a tiff file');
end;
if id(1)~='I'
    error('Not a little-ended tiff file; not supported');
end;

ifdptr=fread(f,1,'uint32');
fseek(f,ifdptr,'bof');
ifdlen=fread(f,1,'uint16');

% Read each of the IFD entries
for i=1:ifdlen
    fseek(f,ifdptr+2+(i-1)*12,'bof');  % go to start of a directory block
    tag=fread(f,1,'*uint16');
    type=fread(f,1,'*uint16');
    count=fread(f,1,'*uint32');
    ptr=fread(f,1,'*uint32');  % pointer to the directory blocks
%     dat=GetData(f,type,count,ptr)
    if tag==37706
        fseek(f,ptr,'bof');
        version=fread(f,1,'*uint32');
        s.comment=fread(f,80,'char')';
        dat1=fread(f,40,'*uint32');  % Read all the integer params from v1
        if version==2  % Hopefully this is what we have
            fseek(f,ptr+584,'bof');
            dv=fread(f,4,'uint8');  % get the date
            time=fread(f,1,'uint32'); % time in seconds since midnight
            dvec=zeros(1,6);        % the Matlab date vector goes here
            dvec(2:3)=dv(2:-1:1);
            dvec(1)=dv(3)+256*dv(4);
            s.datenum=datenum(dvec)+time/(24*3600);  % Store a Matlab datenumber
            s.datestring=datestr(s.datenum);                % Store a string, e.g.
                                                    % 27-Oct-2010 16:45:27

            fseek(f,ptr+2968,'bof');
            pixnm=fread(f,1,'single');
            s.pixA=pixnm*10;        % Angstroms per pixel

            fseek(f,ptr+3272,'bof');
            s.volts=fread(f,1,'single');

            
            fseek(f,ptr+4124,'bof');
            s.ecounts=fread(f,1,'single');  % counts per electron
            
        end;
    end;
end;

fclose(f);


function dat=GetData(f, type, count, pointer)
dat=0;
switch type
    case 5
        fseek(f,pointer,'bof');
        dat=fread(f,count,'uint32');
end;

        
