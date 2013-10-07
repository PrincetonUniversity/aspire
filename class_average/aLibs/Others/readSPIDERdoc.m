function cols = readSPIDERdoc(filename, varargin)
% readSPIDERdoc(filename, [columns])
%  Imports data from a SPIDER document file
%  filename:  file to read
%  varargin: can be either vector of columns to read or
%            a column specified by a single number. 
%    Column 1 is the 1st SPIDER data column.
%    Use column = 0 to get keys.
%    If no columns are specified, then all columns are 
%    returned (without keys).
%
% Examples:
% m = readSPIDERdoc('doc001.dat', 2);        % get 2nd data column.
% m = readSPIDERdoc('doc001.dat', [1 3 5]);  % get columns 1,3,5.
% m = readSPIDERdoc('doc001.dat');           % get all columns.
%
% Returns cols: a matrix of [nrows x ncolumns] 
%
% access by cols(row,column)
% cols(:,i)  gives the ith column
% cols(j,:)  gives the jth row
% cols(j,i)  gives a single element
%
% version 1.0 (Feb 2009) B. Baxter
% Copyright (C) 2009 Health Research Inc.
% Tested in Matlab 7.7.0 (R2008b, R2008a)


if nargin > 1
    arg = varargin{1,1};
    if ~isnumeric(arg)
        disp('readSPIDERdoc: 2nd argument must be column number')
        disp('or vector of column numbers')
        cols = []
        return
    elseif length(arg) == 1
        columns = [arg];
    else
        columns = arg;
    end
else
    columns = -1;
end

% Import the file
newData1 = importdata(filename);

% if there are no comment lines, newData is numeric
if isnumeric(newData1)
    A = newData1;
else
    A = newData1.data;
end

% read all columns
if columns == -1
    [nrows, ncols] = size(A);
    columns = [1:ncols-2];
end

M = [];
% skip the 2nd column that tells number of SPIDER data columns
for c = columns
    if c == 0
        k = 1;
    else
        k = c + 2;
    end
    a = A(:,k);
    M = [M;a'];
end

cols = M';