function assertfield(S,section,fieldname)
%
% ASSERTFIELD   Verify that a field exists in a workflow struct
%
% assertfield(S,section,fieldname)
%   Check that the field S.(section).(fieldname) exists in the struct S.
%
% Yoel Shkolnisky, August 2015

if ~isfield(S, section)
    error('Section %s does not exist is struct %s.',section,inputname(1));
else
    if ~isfield(S.(section),fieldname)
        error('Field %s does not exist is struct %s.',fieldname,section);
    end
end
end
