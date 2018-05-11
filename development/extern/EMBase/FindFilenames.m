function names=FindFilenames(dirName,rexp,isDir)
% given the directory name, pick up filenames that match a pattern rexp.
% Return the names in a cell array.  The optional argument isDir tells the
% function to return directory names instead of file names.
% The pattern follows the Matlab regular expression syntax, which is
% different from Unix.
% Some examples:
% rexp='RawImage_\d*\.tif'  find all names that contain RawImage_81.tif
% rexp='(?i)\d*(u|o)\d*\.(mrc|dm3)' ignore case, find e.g. 01u3000.MRC

if nargin<3
    isDir=0;
end;
d=dir(dirName);
% dirName
% dir
names={};
j=0;
for i=1:numel(d)
    if d(i).isdir==isDir
        q=regexp(d(i).name,rexp);
        if numel(q)>0
            j=j+1;
            names{j}=d(i).name;
        end;
    end;
end;
