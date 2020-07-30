function subsetdata=subsetSTAR(datablock,field,val)
% SUBSETSTAR  Select a subset of STAR data
%
% subsetdata=subsetSTAR(datablock,field,val)
%   Return STAR data containing only records in which field is equal to the
%   given val
%
% Example
%   subsetdata=subsetSTAR(datablock,'rlnClassNumber',1)
%
% Yoel Shkolnisky, December 2017

% NOTE: If STAR files of RELION 3.1 is used, then the structure of the
% STAR file is assumed to contained one optics group (location 1 in the
% stardata array) and one particles group (location 2 in the stardata
% array).
if numel(datablock)==1 % RELION version < 3.1
    idx=zeros(numel(datablock.data),1);
    
    % Set idx to 1 in all indices which satisfies the criterion
    for k=1:numel(datablock.data)
        if isnumeric(val)
            idx(k)=datablock.data{k}.rlnClassNumber==val;
        elseif ischar(val)
            idx(k)=strcmp(datablock.data{k}.(field),val);
        else
            error('Unsupported type for val');
        end
    end
    subsetdata=datablock;
    subsetdata.data=datablock.data(logical(idx));
else
    idx=zeros(numel(datablock(2).data),1);
    
    % Set idx to 1 in all indices which satisfies the criterion
    for k=1:numel(datablock(2).data)
        if isnumeric(val)
            idx(k)=datablock(2).data{k}.(field)==val;
        elseif ischar(val)
            idx(k)=strcmp(datablock(2).data{k}.(field),val);
        else
            error('Unsupported type for val');
        end
    end
    subsetdata=datablock;
    subsetdata(2).data=datablock(2).data(logical(idx));
end
