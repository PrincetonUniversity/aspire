function header=readSPIDERheader(filename)
% readSPIDERheader returns a vector of header values from
% a SPIDER binary file.
%
% version 1.0 (Feb 2009) B. Baxter
% Copyright (C) 2009 Health Research Inc.
% Tested in Matlab 7.7.0 (R2008b, R2008a)

[fp, errmsg] = fopen(filename, 'r');

if fp == -1
  disp(errmsg)
  header = -1;
  return
end

% try reading as big-endian first
nwords = 256;
[h,count] = fread(fp, nwords, 'float32', 'ieee-be');
fclose(fp);

if count < nwords
    %disp('readSPIDERheader: unable to read header');
    header = -1;
    return
end

% Check validity of header.

header = 1;  % assume header is good until shown otherwise
hints = [1 2 5 12 13 22 23];
iforms = [1 3 -11 -12 -21 -22];

% a dummy loop that we can break out of to avoid further tests
for dummy = 1:1 
    % Certain values in the header must be integers
    for i = 1:7
        index = hints(i);
        v = h(index);
        if (v - round(v)) ~= 0
            header = -1;
            break
        end
    end
    
    if header == -1
        break
    end

    % 5th element must be in iforms
    if sum(ismember(iforms,h(5))) == 0
        header = -1;
        break
    end

    % check other values
    labrec = h(13);
    labbyt = h(22);  % total number of bytes in the header
    lenbyt = h(23);
    if labbyt ~= (labrec*lenbyt)
        header = -1;
    end 

end

if header ~= -1
    %s = 'big-endian'
    h(29) = 1;     % unused, put 'endianness' flag here
    header = h;
    return
end 

%====================================================
% try reading as little-endian
header = 1;
fp=fopen(filename, 'r');
[h, count] = fread(fp, nwords, 'float32', 'ieee-le'); 
fclose(fp);

for dummy = 1:1     
    % Certain values in the header must be integers
    for i = 1:7
        index = hints(i);
        v = h(index);
        if (v - round(v)) ~= 0
            header = -1;
            break
        end
    end
    
    if header == -1
        break
    end

    % 5th element must be in iforms
    if sum(ismember(iforms,h(5))) == 0
        header = -1;
        break
    end

    % check other values
    labrec = h(13);
    labbyt = h(22);
    lenbyt = h(23);
    if labbyt ~= (labrec*lenbyt)
        header = -1;
    end 

end

if header ~= -1
    %s = 'little-endian'
    h(29) = 0;
    header = h;
    return
end    

