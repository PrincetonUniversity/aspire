% Test the GPU kernel kernel1.cu
%
% Yoel Shkolnisky, May 2022.

k = parallel.gpu.CUDAKernel('kernel1.ptx','kernel1.cu');

M=10; N=5;
comp_arr=single(randn(M,N)+1i*randn(M,N));
g_comp_arr=gpuArray(comp_arr);
g_real_arr=gpuArray(zeros(M,N,'single'));
g_imag_arr=gpuArray(zeros(M,N,'single'));

blockSize=256;
k.ThreadBlockSize = [blockSize, 1, 1];
k.GridSize = [ceil(M/blockSize), N];

[g_real_arr, g_imag_arr]=feval(k,g_comp_arr,M,N,g_real_arr,g_imag_arr);
real_arr = gather(g_real_arr);
imag_arr = gather(g_imag_arr);

% The following should be precisely zero.
real(comp_arr)-real_arr
imag(comp_arr)-imag_arr