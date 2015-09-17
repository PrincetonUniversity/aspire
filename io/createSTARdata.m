function datablock=createSTARdata(Nrecords,varargin)
% CREATESTARDATA    Create STAR datablock
%
% datablock=createSTARdata(Nrecords,fieldnames)
%   Create an empty STAR data block with the given fieldnames. Records can
%   be added to the datablock using addrecordtoSTARdata. Memory for
%   Nrecords data records is preallcoated (to speed up modifying the
%   structure).
%
% Example:
%   db=createSTARdata('defocus','pixA');
%
% Yoel Shkolnisky, September 2015.

fieldnames=varargin;
datablock.name='data_';
% Create labels
datablock.labels=cell(numel(fieldnames),1);
for k=1:numel(fieldnames)
    datablock.labels{k}=fieldnames{k};
end

% Create data part
if Nrecords<0
    Nrecords=0;
end
datablock.data=cell(Nrecords,1);

% % Create dummy strcuture for initaliztion.
% stripcomment =@(str) regexprep(str,'(\s)*#(.)*',''); % Remove STAR comment 
%         % from a string (comment is in the for '...# comment');
% trimmedlabels=arrayfun(stripcomment,datablock.labels);      
% s=struct;
% for k=1:numel(fieldnames);
%     s.(trimmedlabels{k})=0;
% end
% 
% % Initialize all records.
% for k=1:Nrecords
%     datablock.data{k}=s;
% end
    