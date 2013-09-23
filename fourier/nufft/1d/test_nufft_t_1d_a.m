function test_nufft_t_1d_a
%
% Test the function nufft_t_1d.
%
% Generates random vectors and sampling points of various sizes and
% compares the slow nufft to the fast one.
%
% Prints OK if the error is small, and the ratio between the time of the
% slow and the fast algorithms. Prints ERROR if the error is large.
% 
% Replace 'single' with 'double' below to use double precision accuray.
%
% Yoel Shkolnisky, December 2009.
%

m=128;
n_tests=50;

for k=1:n_tests
    n=min(2^k,2048);
    fprintf('Test for size n=%4d...',n)
    x=rand(m,1); % sampling points
    beta=rand(n,1);
    
    tic
    v1=slow_nufft_t_1d(beta,x);
    t1=toc;
        
    tic;
    v2=nufft_t_1d(beta,x,'single');
    t2=toc;
    
    d=norm(v1-v2)/norm(v1);

    if t1<1.e-12
        t1=0;
        t2=1;
    end

    if d<1.50e-5
        fprintf('OK     d=%7.4e    r=%7.4f\n',d,t1/t2);
    else
        fprintf('ERROR     d=%7.4e    v1=%7.4e     v2=%7.4e\n',d,norm(v1),norm(v2));
    end
end