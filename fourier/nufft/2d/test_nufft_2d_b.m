% Test nufft_2d
%
% Compare nufft_2d to direct computation of the non-equally spaced FFT.
% Sample frequencies are -n/2:n/2-1. alpha is complex.
%
% Yoel Shkolnisky, January 2010.

for n=2:100
    alpha=rand(n)+1i*rand(n);
    [omega_r omega_c]=ndgrid(-n/2:n/2-1,-n/2:n/2-1);
    omega=[omega_r(:) omega_c(:) ];
    % The fast changing frequnecy is omega_r(:) (along the rows), matching
    % the fast dimension of alpha.

    M=2*n;
    
    tic
    f1=slow_nufft_2d(alpha,omega,M);
    t1=toc;

    tic;
    f2=nufft_2d(alpha,omega,'single',M);
    t2=toc;

    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2);
end
