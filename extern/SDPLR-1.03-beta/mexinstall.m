%% User settings

WINLIB1 = [matlabroot '\extern\lib\win32\lcc\libmwlapack.lib'];
WINLIB2 = [matlabroot '\extern\lib\win32\lcc\libmwblas.lib'];

%% Code

% if strcmp(computer,'PCWIN64') | strcmp(computer,'GLNXA64')
%   computer_model = 64; 
% else
%   computer_model = 32; 
% end

matlabversion = sscanf(version,'%f');
matlabversion = matlabversion(1);
matlabversion_1 = floor(matlabversion);
matlabversion_2 = matlabversion - matlabversion_1;
matlabversion_2 = sscanf(num2str(matlabversion_2), '0.%d\n');

% fsp = filesep;

% subdir = strcat('source',fsp);
subdir = 'source/';

fname = cell(15,1);
fname{ 1} = 'mexsdplr.c';
fname{ 2} = 'copystructures.c';
fname{ 3} = 'dataoper.c';
fname{ 4} = 'eigen.c';
fname{ 5} = 'initialize.c';
fname{ 6} = 'lbfgs.c';
fname{ 7} = 'linesearch.c';
fname{ 8} = 'main.c';
fname{ 9} = 'misc.c';
fname{10} = 'params.c';
fname{11} = 'rankreduce.c';
fname{12} = 'readdata.c';
fname{13} = 'sdplrlib.c';
fname{14} = 'timefuncs.c';
fname{15} = 'util.c';

% mexcmd = 'mex -O -v CFLAGS="\$CFLAGS -std=iso9899:1999" -D__MEX ';
mexcmd = 'mex -O CFLAGS="\$CFLAGS -std=iso9899:1999" -D__MEX -largeArrayDims ';

if ispc
 mexcmd = [mexcmd '-D__WIN32 '];
end
if matlabversion_1 >= 7 & matlabversion_2 >= 3
 mexcmd = [mexcmd '-largeArrayDims '];
end
for k = 1:length(fname)
 mexcmd = [mexcmd subdir fname{k} ' '];
end 
mexcmd = [mexcmd 'gsl-1.5/poly/eval.c gsl-1.5/poly/solve_cubic.c '];
if ispc
 mexcmd = [mexcmd '"' WINLIB1 '"' ' ' '"' WINLIB2 '"' ];
else
 mexcmd = [mexcmd '-lmwlapack -lmwblas'];
end

eval(mexcmd)
