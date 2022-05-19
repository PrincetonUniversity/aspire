% Test the kernel2.cu
%
% Yoel Shkolnisky, May 2022.
k = parallel.gpu.CUDAKernel('kernel2.ptx','kernel2.cu');

M=10000; N=10000;
pf1=single(randn(M,N)+1i*randn(M,N));
pf2=single(randn(M,N)+1i*randn(M,N));

g_pf1=gpuArray(pf1);
g_pf2=gpuArray(pf2);

g_corrs=gpuArray(complex(zeros(N,N,'single')));


blockSize=1024;
%k.ThreadBlockSize = [blockSize, blockSize, 1];
%k.GridSize = [ceil(N/blockSize), ceil(N/blockSize)];
k.ThreadBlockSize = [blockSize, 1, 1];
k.GridSize = [ceil(N/blockSize), 1];

disp('Kernel evaluation');
tic;
g_corrs=feval(k,g_pf1,g_pf2,M,N,g_corrs);
corrs_kernel = gather(g_corrs);
t_kernel=toc;

% The following should be of the order of 1.0e-7.
disp('CPU evaluation');
tic;
corrs_ref = pf1.'*conj(pf2);
t_cpu = toc;

% Compare to GPU BLAS
disp('GPU evaluation');
tic;
g_corrs_gpu = g_pf1.'*conj(g_pf2);
corrs_gpu = gather(g_corrs_gpu);
t_gpu = toc;


err=corrs_ref-corrs_kernel;
norm(err(:))/norm(corrs_ref(:))
fprintf("GPU kernel timing = %d secs\n",t_kernel);
fprintf("CPU timing = %d secs\n",t_cpu);
fprintf("GPU BLAS = %d secs\n",t_gpu);
fprintf("kernel speedup = %4.2f\n",t_cpu/t_kernel);