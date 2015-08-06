function tf=fieldexist(S,section,fieldname)
%
% fieldexist   Does a field exists in a workflow struct
%
% fieldexist(S,section,fieldname)
%   Check that the field S.(section).(fieldname) exists in the struct S.
%   return 1 if yes and 0 otherwise.
%
% Yoel Shkolnisky, August 2015

tf=1;
if ~isfield(S, section)
    tf=0;
else
    if ~isfield(S.(section),fieldname)
        tf=0;
    end
end