% Test nufft_1d
%
% Compare nufft_1d to direct computation of the non-equally spaced FFT.
% Vectors lengths are powers of two. Frequencies are -n/2:n/2-1.
%
% Yoel Shkolnisky, January 2010.

for k=2:11
    n=2^k;
    alpha=rand(n,1);
    omega=(-n/2:n/2-1).';
    
    tic
    f1=slow_nufft_1d(alpha,omega);
    t1=toc;
        
    tic;
    f2=nufft_1d(alpha,omega,'single');
    t2=toc;
    
    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2);
end
