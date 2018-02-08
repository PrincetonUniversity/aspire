function kv=struct2keyvalue(s)
%
% STRUCT2KEYVALUE   Convert a strucut to a cell array of key-value
%
% kv=strucut2keyvalue(s)
%   Convert the struct s to a cell array of key/value pairs, where in each
%   pair the first item is the field name and the second item is the value
%   of the field.
%
% Example:
%  s.f1=1; s.f2='a';
%  kv=struct2keyvalue(s)
%
% Yoel Shkolnisky, June 2017

fields=fieldnames(s);
kv=cell(numel(fields)*2,1);
for i=1:numel(fields)
    key=fields{i};
    val=s.(key);
    kv{2*i-1}=key;
    kv{2*i}=val;
end
