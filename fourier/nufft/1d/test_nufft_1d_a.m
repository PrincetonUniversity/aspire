% Test nufft_1d
%
% Compare nufft_1d to the ifft.
%
% Yoel Shkolnisky, January 2010.

for k=2:15
    n=2^k;
    alpha=rand(n,1);
    omega=(-n/2:n/2-1).';
    
    tic
    f1=nufft_1d(alpha,omega,'single');
    t1=toc;
    
    tic;
    f2=n*icfft(alpha);
    t2=toc;
    
    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2);
end
