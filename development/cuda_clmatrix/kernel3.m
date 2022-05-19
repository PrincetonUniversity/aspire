% Test kernel3.cu
%
% Yoel Shkolnisky, May 2022.

k = parallel.gpu.CUDAKernel('kernel3.ptx','kernel3.cu');

N=single(5000);
rots=single(rand_rots(N));
n_theta=single(360);

% Generate reference clmatrix
disp('CPU evaluation');
tic;
clmatrix_ref=clmatrix_cheat(rots,n_theta);
t_cpu=toc;

% Prepare output array
disp('Kernel evaluation');
tic;
g_rots = gpuArray(single(rots));
g_clmatrix=gpuArray(zeros(N,N,'single'));

% Call GPU kernel
blockSize=1024;
k.ThreadBlockSize = [blockSize, 1, 1];
k.GridSize = [ceil(N/blockSize), N];

g_clmatrix=feval(k,g_rots,N,n_theta,g_clmatrix);
clmatrix=gather(g_clmatrix);
t_kernel=toc; 

% Test the results.
% Allow of a difference of up to 1 line, due to different roundoff.
accuracy = comparecl(clmatrix,clmatrix_ref,n_theta,360/n_theta*1.5); % Allow error of up to one and a half lines
if accuracy<1
    error('Errors too large. Accuracy = %3.2f',accuracy);
end

fprintf("GPU kernel timing = %d secs\n",t_kernel); % This timing includes data transfer to the GPU.
    % Note that data transfer takes the overwhelming majority of time.
fprintf("CPU timing = %d secs\n",t_cpu);
fprintf("kernel speedup = %4.2f\n",t_cpu/t_kernel); 
