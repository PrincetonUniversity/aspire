mex -O nufft3dauxmx.cpp -I../common
mex -O nufft3dauxmx_ref.cpp -I../common
mex nufftt3dexecutemx_1.cpp ../common/printarr.cpp -I../common/
mex nufftt3dexecutemx_2.cpp ../common/printarr.cpp -I../common/