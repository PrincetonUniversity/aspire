function pathstr=AddSlash(pathstr)
% function pathstr=AddSlash(pathstr)
% given a path string, check that it ends with a slash.  If not, add one.

n=numel(pathstr);
if n>0 && pathstr(n) ~='/'
    pathstr=[pathstr '/'];
end;
