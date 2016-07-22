function datablock=addrecordtoSTARdata(datablock,idx,varargin)
% ADDRECORDTOSTARDATA   Add a new record to a STAR data block
%
% datablock=addrecordtoSTARdata(datablock,idx,varargin)
%   Add the given input as a new record to the input datablock at index
%   idx. The modified datablock is returned. The input is in the for of
%   'key','value','key','value',... containing all the required fields and
%   their values. All fields specified by the field labels of datablock
%   must be provided. if idx is -1, the new record is added at the end.
%
% Example:
%   datablock=addrecordtoSTARdata(datablock,'voltage',300,...
%       'DefocusU',2.3344699219e+03,'DefocusV',2.3445949219e+03',...
%       'DefocusAngle',36.7,'Cs',2.0,'pixA',1.4,'A',0.1);
%
% Yoel Shkolnisky, September 2015.

% Strip comments from labels
stripcomment =@(str) regexprep(str,'(\s)*#(.)*',''); % Remove STAR comment 
        % from a string (comment is in the for '...# comment');
trimmedlabels=arrayfun(stripcomment,datablock.labels);      

% Verify that all fields were provided
fieldnames=varargin(1:2:end);
vals=varargin(2:2:end);

if numel(fieldnames)~=numel(vals)
    error('Incompatible number of keys and values.');
end

if ~isempty(setxor(trimmedlabels,fieldnames))
    labelsdiff=setdiff(trimmedlabels,fieldnames);
    error('Not all required fields were provided. Missing fields:\n%s',...
        sprintf('%s\n', labelsdiff{:}));
end

% Update data
if idx==-1
    ndata=numel(datablock.data);    
    idx=ndata+1;
end

%datablock.data{idx}=struct(varargin{:});

for k=1:numel(fieldnames)
    datablock.data{idx}.(fieldnames{k})=vals{k};
end

