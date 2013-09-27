function b=StructsEqual(A,B)
% function b=StructsEqual(A,B)
% Test whether structs A and B have the same field names and values.  Here
% fields are assumed to be scalars or arrays, but not structs.
% Someday we should overload the 'eq' function.
fa=sort(fieldnames(A));
fb=sort(fieldnames(B));
na=numel(fa);
b=0;
if na ~= numel(fb)
    return
end;
for i=1:na
    fname=char(fa(i));
    if ~(strcmp(fname,char(fb(i))) && ...
        all(size(A.(fname))==size(B.(fname))) && ...
        all(A.(fname)(:)==B.(fname)(:)));
        return
    end;
end;
b=1;
