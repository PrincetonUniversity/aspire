function fnameext=addExtIfNeeded(fname,ext)
%
% ADDEXTIFNEEDED    Add extension to filename if needed.
%
% fname=addExtIfNeeded(fname,ext)
%   If filename fname does not haven the extension ext then appened ext to
%   fname. Else return fname. Extension must be prefixed by ".".
%
% Examples:
%       addExtIfNeeded('vol.mrc','.mrc')
%   returns 'vol.mrc'.
%       addExtIfNeeded('vol','.mrc')
%   returns 'vol.mrc'
%
% Yoel Shkolnisky, June 2017.

if nargin~=2
    error('filename and extension are required');
end

[d,n,e]=fileparts(fname);
if isempty(e)
    fnameext=fullfile(d,[n ext]);
else
    fnameext=fname;
end
    
    
