function mprintf(handles,varargin)
% function mprintf(handles,varargin)
% Write the same formatted text to multiple files.  handles is an array of
% file handles.  The function fprintf is called multiple times, once for
% each handle number.  Handles of 0 are ignored.
for i=1:numel(handles)
    if handles(i)>0
        fprintf(handles(i),varargin{:});
    end;
end
