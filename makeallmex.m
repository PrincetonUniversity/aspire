function makeallmex
% MAKEALLMEX   Compile all ASPIRE mex files
%
%   MAKEALLMEX
%
% The function should be called only once upon installation of ASPIRE to
% compile all mex files for the traget platform.
%
% Yoel Shkolnisky, September 2013.

currPath=pwd;
[aspirePath,~,~]=fileparts(mfilename('fullpath'));
traversedirrectorytree(aspirePath,@runmakemex);
cd(currPath);

function runmakemex(fullFileName)
[mexPath,fname]=fileparts(fullFileName);
if strcmp(fname,'makemex')
    cd(mexPath);
    fprintf('Running %s\n',fullFileName);
    run(fname);
end
    
    

function traversedirrectorytree(path,fileOperation)
%
% TRAVERSEDIRRECTORYTREE Apply "fileOperation" to a directory tree.
%
%   TRAVERSEDIRRECTORYTREE(path,operation) traverse a directory tree whose
%   root is path and apply "fileOperation" to each node. The function handle
%   operation has the signature operation(fullFileName).
%
% Example:
%  dispfile= @(str) disp(str);
%  traversedirrectorytree(path,dispfile) 
%
% Yoel Shkolnisky, September 2013

if ~strcmp(path(end), filesep)
    path(end+1)=filesep;
end
dirInfo= dir(path);
files=~[dirInfo.isdir];
fileNames={dirInfo(files).name};
%disp(path);
if ~isempty(fileNames)
    for i=1:length(fileNames)
        % Do whathever (show file names?)
        fileOperation(char(fullfile(path,fileNames(i))));
    end
end

% For each subdir, call itself again
isDir=[dirInfo.isdir];
dirNames={dirInfo(isDir).name};
dirNames(strcmp(dirNames, '.') | strcmp(dirNames, '..'))=[];

for i=1:length(dirNames)
    traversedirrectorytree([path dirNames{i} filesep],fileOperation);    
end