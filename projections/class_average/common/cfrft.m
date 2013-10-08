% function w=cfrft(x,alpha)
%
% Aliased fractional Fourier transform of the sequence x.
% The FRFT is computed using O(nlogn) operations.
%
% x       The sequence whose FRFT should be computed. Can be of odd or even
%         length. Must be a 1-D row vector.
% alpha   The parameter alpha of the fractional Fourier transform.
%
% Returns the aliased FRFT with parameter alpha of the sequence x.
% The fractional Fourier transform w of the sequence x (with parameter alpha) is defined by
%                   n/2-1
%       w(k) =       sum  x(u)*exp(-2*pi*i*k*u*alpha/N),  -n/2 <= k <= n/2-1, N=length(x).
%                   u=-n/2
%
% 
% This function is the same as cfrftV2. It uses the less padding (3m as in the paper)
% and therefore it is more memory efficient. This may cause the lengths of the sequences to be non-optimal 
% for Matlab's FFT function. The function cfrftV2 uses more memory (for padding) but uses FFTs
% of dyadic length.
%
% Yoel Shkolnisky 22/10/01


function w=cfrft(x,alpha)
m=length(x);
j=lowIdx(m):hiIdx(m);
j2= lowIdx(2*m):hiIdx(2*m);
E=i*pi*alpha;

y=x.*exp(-E*j.^2/m);
y=[zeros(1,m),y,zeros(1,m)];

z=zeros(1,3*m);
l=toUnaliasedIdx(-m,3*m);
z(l:l+length(j2)-1)=exp(E*j2.^2/m);

Y=cfft(y);
Z=cfft(z);
W=Y.*Z;
w=icfft(W);
w=w(toUnaliasedIdx(lowIdx(m),3*m):toUnaliasedIdx(hiIdx(m),3*m));
w=w.*exp(-E*j.^2/m);
