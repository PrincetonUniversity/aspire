% Some tests of kernel4.cu
% Yoel Shkolnisky, May 2022.

% k = parallel.gpu.CUDAKernel('kernel4.ptx','kernel4.cu');
% 
% N=360; M=128; P=129;
% A=single(randn(N,P)+1i*randn(N,P));
% B=single(randn(P,M)+1i*randn(P,M));
% 
% g_A=gpuArray(complex(A));
% g_B=gpuArray(complex(B));
% g_C=gpuArray(complex(zeros(N,M,'single')));
% 
% 
% blockSize=1;
% k.ThreadBlockSize = [1, 1, 1];
% k.GridSize = [1, 1];
% 
% disp('Kernel evaluation');
% tic;
% g_C=feval(k,g_A,g_B,g_C,N,M,P);
% t_kernel=toc;
% C = gather(g_C);
% 
% 
% tic;
% Cref=A*B;
% t_cpu = toc;
% 
% err=C-Cref;
% disp(norm(err(:))/norm(C(:)))
% 
% fprintf("GPU kernel timing = %d secs\n",t_kernel); % This timing includes data transfer to the GPU.
%     % Note that data transfer takes the overwhelming majority of time.
% fprintf("CPU timing = %d secs\n",t_cpu);
% fprintf("kernel speedup = %4.2f\n",t_cpu/t_kernel); 
% 
% k = parallel.gpu.CUDAKernel('kernel4.ptx','kernel4.cu');
% n=300;  m=500;
% A=single(randn(n,m)+1i*randn(n,m));
% V=single(randn(n,1)+1i*randn(n,1));
% 
% g_A=gpuArray(A);
% g_V=gpuArray(V);
% g_AV=gpuArray(complex(zeros(n,m,'single')));
% 
% 
% blockSize=1;
% k.ThreadBlockSize = [1, 1, 1];
% k.GridSize = [1, 1];
% 
% g_AV=feval(k,g_A,g_V,g_AV,n,m);
% t_kernel=toc;
% AV = gather(g_AV);
% 
% AVref=bsxfun(@times,A,V);
% norm(AV-AVref,'fro')/norm(AV,'fro')


k = parallel.gpu.CUDAKernel('kernel4.ptx','kernel4.cu');
n=5000;  m=1000;
A=single(randn(n,m));

g_A=gpuArray(single(A));
g_val=gpuArray(single(zeros(1)));
g_idx=gpuArray(single(zeros(1)));


blockSize=1;
k.ThreadBlockSize = [1, 1, 1];
k.GridSize = [1, 1];

tic;
[g_val,g_idx]=feval(k,g_A,m*n,g_val,g_idx);
val=gather(g_val);
idx=gather(g_idx);
t_kernel=toc;

tic;
[vv,ii]=max(A(:));
t_cpu=toc;

assert(vv==val);
assert(idx==ii-1);

fprintf("GPU kernel timing = %d secs\n",t_kernel); % This timing includes data transfer to the GPU.
fprintf("CPU timing = %d secs\n",t_cpu);
fprintf("kernel speedup = %4.2f\n",t_cpu/t_kernel); 
