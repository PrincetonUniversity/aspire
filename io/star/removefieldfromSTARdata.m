function datablock=removefieldfromSTARdata(datablock,fieldname)
% REMOVEFIELDFROMSTARDATA    Remove field from a STAR structure
%
% datablock=removefieldfromSTARdata(datablock,fieldname)
%   Remove the given fieldname from the STAR data srtucture.The field is
%   remove both from tjhe labels section as well as from the data section.
%   Returns the modified data block.
%   
% Yoel Shkolnisky, September 2015.

% Remove the field from the labels section.
for k=numel(datablock.labels):-1:1
    idx=regexp(datablock.labels{k},[fieldname,'(\s)*#(.)*'],'once');
    if ~isempty(idx)  
        datablock.labels(k)=[];
    end;
end

% Remove the field from the data section.
ndata=numel(datablock.data);
for k=1:ndata
    datablock.data{k}=rmfield(datablock.data{k},fieldname);
end
