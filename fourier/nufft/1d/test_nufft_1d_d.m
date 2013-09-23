% Test nufft_1d
%
% Similar to test_nufft_1d_c.m but the number of output points is
% different from the number of input points. The number of output points is
% random.
%
% We sometimes get errors that are slightly bigger (3.0e-5) only if N is
% very small and M>>N.
% 
% Yoel Shkolnisky, January 2010.

rand('state',111);
K=1:1000;
N=ceil(rand(numel(K),1)*1024);
N=sort(N);
M=ceil(rand(numel(K),1)*1024);

for k=K;
    n=N(k);
    alpha=rand(n,1);
    omega=(rand(n,1)-1/2)*n;
    
    tic
    f1=slow_nufft_1d(alpha,omega,M(k));
    t1=toc;
        
    tic;
    f2=nufft_1d(alpha,omega,'single',M(k));
    t2=toc;
    
    checkerror(n,f1,f2,alpha,1.0e-5,t1,t2);
end
