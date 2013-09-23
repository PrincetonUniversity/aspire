% Test nufft_2d
%
% Compare nufft_2d to the ifft2.
%
% Yoel Shkolnisky, January 2010.

for k=2:8 
    n=2^k;
    alpha=rand(n)+1i*rand(n);
    [omega_r omega_c ]=ndgrid(-n/2:n/2-1,-n/2:n/2-1);
    omega=[omega_r(:) omega_c(:) ];
    % The fast changing frequnecy is omega_r(:) (along the rows), matching
    % the fast dimension of alpha.

    tic
    f1=nufft_2d(alpha,omega,'single',n);
    t1=toc;

    tic;
    f2=n^2*icfft2(alpha);
    t2=toc;

    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2)
       
end
