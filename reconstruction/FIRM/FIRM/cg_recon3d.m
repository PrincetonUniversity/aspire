function [x, err, iter, flag] = cg_recon3d(A, b, tol, max_it, x0)
% Fastest version.
%
% Input:
%   A       convolution kernel
%   b       backprojection of the images
%   tol     error toleration
%   max_it  maximun number of iterations
%   x0      initial guess
% Output:
%   x       solution
%   err     estimated error
%   iter    number of iterations
%   flag    INTEGER: 0 = solution found to tolerance
%                    1 = no convergence given max_it
%
% Lanhui Wang, Princeton University, Oct 5, 2010
n=size(b,1);

if isempty(tol), tol=1e-3; end
if isempty(max_it), max_it=n^3; end
if isempty(x0), x0 = zeros(n,n,n); end



%%
% initialization
x=x0;
flag = 0;
iter = 0;
bnrm2 = norm( b(:) );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

% NEW definition for the SAME varibles
% a (of size 2n) is now the Fourier transfrom of the first column
% of the circulant matrix [A X;X A]
a=zeros(2*n,2*n,2*n);
a(1:n,1:n,1:n)=A(n:2*n-1,n:2*n-1,n:2*n-1);
a(n+1:2*n,1:n,1:n)=A([n 1:n-1],n:2*n-1,n:2*n-1);
a(1:n,n+1:2*n,1:n)=A(n:2*n-1,[n 1:n-1],n:2*n-1);
a(n+1:2*n,n+1:2*n,1:n)=A([n 1:n-1],[n 1:n-1],n:2*n-1);
a(1:n,1:n,n+1:2*n)=A(n:2*n-1,n:2*n-1,[n 1:n-1]);
a(n+1:2*n,1:n,n+1:2*n)=A([n 1:n-1],n:2*n-1,[n 1:n-1]);
a(1:n,n+1:2*n,n+1:2*n)=A(n:2*n-1,[n 1:n-1],[n 1:n-1]);
a(n+1:2*n,n+1:2*n,n+1:2*n)=A([n 1:n-1],[n 1:n-1],[n 1:n-1]);
a = fftn(a);



% Compute Toeplitz matrix-vector product by fft
%   w = [x;zeros(n,n,n)];
w=zeros(2*n,2*n,2*n,class(b));
w(1:n,1:n,1:n)=x;
Aw = ( ifftn( a .* fftn(w) ) ); Aw = Aw(1:n,1:n,1:n);
r = b - Aw;

error=norm(r(:));
if ( norm(b(:)) == 0 )
    os = sprintf(['The right hand side vector is all zero so CPCG '...
        'returned an all zero solution without iterating.']);
    disp(os);
    return,end
if ( error < tol )
    os = sprintf(['The initial guess is the solution with relative residual '...
        '%0.2g so CPCG without iteration'],error);
    disp(os);
    return, end

err=zeros(1,max_it);

for iter = 1:max_it                       % begin iteration
    z  = r;
    rho = (r(:)'*z(:));
    
    if ( iter > 1 ),                       % direction vector
        beta = rho / rho_1;
        p = z + beta*p;
    else
        p = z;
    end
    w(1:n,1:n,1:n)=p;
    Aw = ( ifftn( a .* fftn(w) ) );
    Aw=Aw(1:n,1:n,1:n);
    
    q = Aw;
    alpha = rho / (p(:)'*q(:) );
    x = x + alpha * p;                    % update approximation vector
    
    r = r - alpha*q;                      % compute residual
    error=norm(r(:));               % check convergence
    
    
    
    err(iter)=error/bnrm2;
    if ( error/bnrm2 <= tol ),
        os = sprintf(['CPCG converged at iteration %d to a solution with '...
            'relative residual %0.2g'], iter,error/bnrm2);
        disp(os);
        break, end
    
    rho_1 = rho;
    if mod(iter,10)==1
        fprintf('iter %d  vol error:%f\n',iter,err(iter));
    end
    
    
end

if ( error/bnrm2 > tol ),
    os = sprintf(['CPCG stopped after %d iteration'...
        'and has relative residual %0.2g'],max_it,error/bnrm2);
    flag = 1; disp(os);                   % no convergence
end

