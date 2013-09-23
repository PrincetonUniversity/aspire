% Test nufft_2d
%
% Compare nufft_2d to direct computation of the non-equally spaced FFT.
% Sample frequencies are random (see test_nufft_2d_b.m).
%
% Yoel Shkolnisky, January 2010.

rand('state',111);
for n=2:100
    alpha=rand(n)+1i*rand(n);
    omega=(rand(n^2,2)-1/2)*n;

    M=n;

    tic
    f1=slow_nufft_2d(alpha,omega,M);
    t1=toc;

    tic;
    f2=nufft_2d(alpha,omega,'single',M);
    t2=toc;

    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2);
end
