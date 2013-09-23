% function y=cfft2e(x)
%
% Aliased FFT2 of the 2D array x, assuming that both the spatial sampling
% points and the frequnecy sampling points are symmetric around 0.
% The FFT is computed using O(nlogn) operations.
%
% x   The 2D sequence whose FFT should be computed. 
%     Can be of odd or even length. Must be a 2D array.
%
% Returns the aliased FFT of the sequence x.
%
% Yoel Shkolnisky, September 2013.

function y=cfft2e(x)

[m,n]=size(x);

if mod(m,2)==1 || mod(n,2)==1;
    error('Input must have even dimensions');
end

% Pre-multiply by appropriate weights
k1=-m/2:m/2-1;
k2=-n/2:n/2-1;

f1=exp(-2.*pi*1i.*k1/(2*m)).';
f2=exp(-2.*pi*1i.*k2/(2*n));

x=bsxfun(@times,x,f1);
x=bsxfun(@times,x,f2);

y=cfft2(x);

% Post-multiply
w1=2.*pi.*(k1+1/2)/m;
w2=2.*pi.*(k2+1/2)/n;

g1=exp(-1i.*w1/2).';
g2=exp(-1i.*w2/2);

y=bsxfun(@times,y,g1);
y=bsxfun(@times,y,g2);
