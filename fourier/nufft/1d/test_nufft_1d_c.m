% Test nufft_1d
%
% Compare nufft_1d to direct computation of the non-equally spaced FFT.
% Vectors have random lengths. Frequencies are random.
%
% Yoel Shkolnisky, January 2010.

rand('state',111);
K=1:100;
N=ceil(rand(numel(K),1)*1024);
N=sort(N);
for k=K;
    n=N(k);
    alpha=rand(n,1);
    omega=(rand(n,1)-1/2)*n;
    
    tic
    f1=slow_nufft_1d(alpha,omega);
    t1=toc;
        
    tic;
    f2=nufft_1d(alpha,omega,'single');
    t2=toc;
    
    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2);
end
