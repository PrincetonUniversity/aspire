% Test nufft_3d
%
% Compare nufft_3d to the ifft3.
%
% Yoel Shkolnisky, January 2010.

n=32; % Must be even in this test.
alpha=rand(n,n,n)+1i*rand(n,n,n);
[omega_x omega_y omega_z]=ndgrid(-n/2:n/2-1,-n/2:n/2-1,-n/2:n/2-1);
omega=[omega_x(:) omega_y(:)  omega_z(:)];

tic
f1=nufft_3d(alpha,omega,'single',n);
t1=toc;

tic;
f2=n^3*icfftn(alpha);
t2=toc;

checkerror(n,f1,f2,alpha,1.0e-5,t1,t2)