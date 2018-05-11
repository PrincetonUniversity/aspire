function datablock=addfieldtoSTARdata(datablock,fieldname,fieldvals)
% ADDFIELDTOSTARDATA    Add new field to a STAR structure
%
% datablock=addfieldtoSTARdata(datablock,fieldname,fieldvals)
%   Add the given fieldname to the STAR data srtucture. If fieldvals is a
%   scalar, then all records get the same value for the new field.
%   Otherwise, length of fieldvals must be equal to the number of data
%   records in the datablock struct. Currently fieldvals must be a real
%   scalar or vector. Returns the modified data block.
%
% Example:
% CTFdata=readSTAR('ctffile.star');
% CTFdata=addfieldtoSTARdata(CTFdata,'pixA',1.23)
% writeSTAR(CTFdata,'ctffile_new.star');
%   
% Yoel Shkolnisky, August 2015.

% Add the new field to the list of labels.
nlabels=numel(datablock.labels);
datablock.labels{nlabels+1}=fieldname;

% Update all records in the data section.
ndata=numel(datablock.data);

if ~isreal(fieldvals)
    error('Values to add must be real. No complex numbers are currently allow.');
end
if numel(fieldvals)==1 % We have a scalar
    fieldvals=repmat(fieldvals,ndata,1);    
end

if size(fieldvals,1)~=ndata
    error('Length of values must be either 1 or equal to number of data records.');
end

for k=1:ndata
    datablock.data{k}.(fieldname)=fieldvals(k);
end

