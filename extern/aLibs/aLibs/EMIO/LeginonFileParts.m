function [path basename imtype ext]=LeginonFileParts(name)
% function [path basename imtype ext]=LeginonFileParts(name)
% Similar to fileparts, but strips off final alpha characters and returns
% them as imtype.  Unlike fileparts, the '/' in the path is kept.
% For example, suppose name='path/10sep19a_a_00006gr_00017sq_v01_00021hl_v02_00004em.mrc'
% Then the returned values will be
%   path='path/'
%   basename='10sep19a_a_00006gr_00017sq_v01_00021hl_v02_00004'
%   imtype='em'
%   ext='.mrc'
% and the original name is obtained simply by concatenation.

[path name1 ext]=fileparts(name);
if numel(path)>0
    path=[path filesep];
end;
len1=numel(name1);

imtype='';
% Search through all alpha characters at the end of the name
i=len1;
while i>0 && name1(i)>='A' && name1(i)<='z'
    imtype=[name1(i) imtype];
    i=i-1;
end;

if i<1  % No non-alpha characters found; return the whole name as basename.
    imtype='';
    basename=name1;
else
    basename=name1(1,1:i);
end;
