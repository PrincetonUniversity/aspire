% /usr/bin/scl enable devtoolset-2 mex -O nufft2doptimizedmx.cpp -I../common
% /usr/bin/scl enable devtoolset-2 mex -O nufft2doptimizedmx.cpp -I../common
% /usr/bin/scl enable devtoolset-2 mex -O nufft2doptimizedmx.cpp -I../common
% /usr/bin/scl enable devtoolset-2 mex -O nufft2doptimizedmx.cpp -I../common

disp('Enable scl ifusing MATLAB 2015 since gcc 4.8 unavailable')

mex -O nufft2doptimizedmx.cpp -I../common
mex -O nufftt2dexecutemx.cpp -I../common
mex -O nufftt2dnoqmx.cpp -I../common
mex -O nufftt2dpreparemx.cpp -I../common