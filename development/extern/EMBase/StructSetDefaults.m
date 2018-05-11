function s=StructSetDefaults(names,values,s)
% For each field that is absent, create the field and insert the value.
for i=1:numel(names)
    name=names{i};
    if ~isfield(s,name)
        s=setfield(s,name,values{i});
    end;
end;
