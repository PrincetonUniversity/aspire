function test_nufft_t_3d_c
%
% Test the functions nufft_t_3d_prepare_2 and nufft_t_3d_execute_2.
% See test_nufft_t_3d_a for more information.
%
% Yoel Shkolnisky, February 2010.
%

m=1024;
n_tests=50;

for k=1:n_tests
    n=min(2^k,64);
    fprintf('Test for size n=%4d...',n)
    x=rand(m,3); % sampling points
    beta=rand(n,n,n);      
    
    tic
    v1=nufft_t_3d(beta,x);
    t1=toc;
           
    tic;
    prepdata=nufft_t_3d_prepare_2(beta,'single');
    v2=nufft_t_3d_execute_2(x,prepdata);
    t2=toc;
    
    d=norm(v1-v2)/norm(v1);

    if t1<1.e-12
        t1=0;
        t2=1;
    end
       
    if d<1.0e-14
        fprintf('OK    d=%7.4e    r=%7.4f\n',d,t1/t2);
    else
        fprintf('ERROR\n');
    end
end