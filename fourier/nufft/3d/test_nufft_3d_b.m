% Test nufft_3d
%
% Compare nufft_3d to direct computation of the non-equally spaced FFT.
% Frequencies are -n/2:n/2-1.
%
% Yoel Shkolnisky, January 2010.

for n=2:10
    alpha=rand(n,n,n)+1i*rand(n,n,n);
    [omega_x omega_y omega_z]=ndgrid(-n/2:n/2-1,-n/2:n/2-1,-n/2:n/2-1);
    omega=[omega_x(:) omega_y(:)  omega_z(:)];
    % The fast changing frequnecy is omega_r(:) (along the rows), matching
    % the fast dimension of alpha.
    
    M=n;
    
    tic
    f1=slow_nufft_3d(alpha,omega,M);
    t1=toc;
    
    tic;
    f2=nufft_3d(alpha,omega,'single',M);
    t2=toc;
    
    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2)
end
