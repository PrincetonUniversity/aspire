function [fields vals]=GetStructValues(s)
% Get field names and values from structure s.
% Like the built-in function, but returns values as well, both as cell
% arrays.
fields=fieldnames(s);
np=numel(fields);
vals=cell(np,1);
for i=1:np
    vals{i}=s.(fields{i});
end;
