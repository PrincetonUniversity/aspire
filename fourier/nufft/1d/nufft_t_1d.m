function g=nufft_t_1d(beta,x,precision)
%
% Compute the adjoint non-equally spaced FFT.
%
% The function computes the sums
%               n/2
%        g(j) = sum beta(k)*exp(i*k*x(j))
%              k=-n/2
% for j=1,...,n.
%
% The complexity of the algorithm is O(n*log(1/eps)+m*n(log(n)))
%
% m is the oversampling factor (m=2).
%
% Input parameters:
%    beta     Coefficients in the sums above. Real or complex numbers.
%    x        Sampling points. Real numbers in the range [-pi,pi]. 
%    precision  'double' or 'single'. Default is 'single'.
%
%    beta and x need not be the same length, but x should be longer than
%    beta.
%
% Output:
%    g        The sums defined above.
%
% Yoel Shkolnisky, December 2006.
%
% Revisions:
% Filename changed from nufft_t to nufft_t_1d. Fixed editor warnings.
% Change default for precision to single. (Y.S. December 22, 2009) 

if nargin<3
    precision='single';
end

if (~strcmpi(precision,'single')) && (~strcmpi(precision,'double'))
    precision='single';
    warning('MATLAB:unkownOption','Unrecognized precsion. Using ''single''.');
end

if strcmpi(precision,'double')
    b=1.5629;
    m=2;
    q=28;
else %single precision
    b=0.5993;
    m=2;
    q=10;
end


if mod(q,2)==1
    error('Set q to an even integer');
end

if length(size(beta))~=2
    error('beta must be a vector');
end

if (size(beta,1)~=1) && (size(beta,2)~=1)
    error('beta must be a vector');
end

if length(size(x))~=2
    error('x must be a vector')
end

if size(x,2)~=1
    error('x must be a vector');
end

n=length(beta);
len_x=length(x);

% if len_x<n
%     error('x must be longer than beta');
% end

l1=ceil((n-1)/2);
h1=floor((n-1)/2);
l2=ceil((len_x-1)/2);
h2=floor((len_x-1)/2);
beta=[zeros(l2-l1,1) ; beta ; zeros(h2-h1,1)];

x=x(:);
beta=beta(:);

n=length(beta);
low_idx=-ceil((n-1)/2);
high_idx=floor((n-1)/2);
idx=low_idx:high_idx;

Q=zeros(len_x,q+1);
%E=zeros(n,1);
g=zeros(len_x,1);

nu=round(x*m*n/(2*pi));

for k=-q/2:q/2
    Q(:,k+q/2+1)=exp(-(x.*m.*n./(2*pi)-(nu+k)).^2/(4*b))/(2*sqrt(b*pi));
end


E=exp(b.*(2.*pi.*idx./(m*n)).^2);
E=E(:);
u=beta.*E;

low_idx_u=ceil((m*n-1)/2);
high_idx_u=floor((n*m-1)/2);
u=[zeros(low_idx_u+low_idx,1) ; u ; zeros(high_idx_u-high_idx,1)];
U=fftshift(ifft(ifftshift(u)));
U=U*length(U);

for k=-q/2:q/2
    idx=nu+k;
    idx=mod(idx+low_idx_u,m*n);
    g=g+Q(:,k+q/2+1).*U(idx+1);
end



