function writeSPIDERdoc(filename, A, varargin)
% Write data out in the format of a SPIDER document file.
% Column headings may be specified.
%
% Examples
% writeSPIDERdoc('doc001.dat', A)
% writeSPIDERdoc('coords.snd', A, {'X axis'; 'Y axis'})
%
% Where A is a 2D matrix of float values.
%
% Optional 3rd arg is a cell array of column header strings,
% created by e.g., hdrs = { 'column1'; 'column2'; 'column3'}
% Note the curly braces and semicolons.
%
% SPIDER files do not have a specific file extension.
%
% version 1.0 (Feb 2009) B. Baxter
% Copyright (C) 2009 Health Research Inc.
% Tested in Matlab 7.7.0 (R2008b, R2008a)

[fp, errmsg] = fopen(filename, 'w');
if fp == -1
  disp(errmsg)
  return
end

hdrstr = makeDocfileHeader(filename);
fprintf(fp,'%s\n', hdrstr);

if nargin > 2
    %celldisp(varargin)
    colheaders = makeColumnHeaders(varargin);
    %disp(colheaders)
    fprintf(fp,'%s\n', colheaders);
end

ndimsa = ndims(A);

if ndimsa > 2
  disp('writeSPIDERdoc: input matrix must be 1 or 2 dimensions')
  fclose(fp);
end

[nrows ncols] = size(A);

for key = 1:nrows
    s = sprintf('%5d %2d', key, ncols);
    for i = 1:ncols
        s1 = sprintf('%11g',A(key,i));
        s = [s s1];
    end
    %disp(s)
    fprintf(fp,'%s\n', s);
end

fclose(fp);

% ===========================================================
% If the array of header strings is created by the command:
%    hdrs = { 'column1'; 'column2'; 'column3'}
% then varargin is a cell array and the results of 
% celldisp(varargin) look something like:
%    varargin{1}{1} = AMPLITUDE
%    varargin{1}{2} = X-AXIS
%    varargin{1}{3} = SPAT.FREQ
%
% The 1st element of varargin (right-hand index) is a cell
% array of strings, where each element is a character array:
% h =  varargin{1,1}
% celldisp(h)
%    h{1} = AMPLITUDE
%    h{2} = X-AXIS
%    h{3} = SPAT.FREQ
% ischar(h{1}) == 1   % the h{i} are your column header strings.
%

function docstr=makeColumnHeaders(aCellArray)

docstr = '';

if iscell(aCellArray) && iscellstr(aCellArray{1,1})
    headers = aCellArray{1,1};
    numheaders = length(headers);
    %disp(sprintf('there are %d headers:', numheaders));
    %for i =1:numheaders
    %   disp(sprintf('  %s', headers{i}));
    %end
else
    return
end

docstr = ' ; /    ';

for i =1:numheaders
    s = sprintf('%11s', headers{i});
    docstr = [docstr s];
end

% ===========================================================
%
%
function hdr=makeDocfileHeader(filename)

[pathstr, basename, dext, ver] = fileparts(filename);

name = [basename dext];

ne = length(dext);
ext = dext(2:ne);  % remove the dot from the extension

idate = date;
t = clock;  % YY MM DD HH MM SS(not an int)
hr  = sprintf('%02d', t(4));
min = sprintf('%02d', t(5));
sec = sprintf('%02d', round(t(6)));
itime = strcat(hr,':',min,':',sec);

hdr = [' ;mat/' ext '   ' idate ' AT ' itime '   ' name];
