% function w=cfrftV2(x,alpha)
%
% Aliased fractional Fourier transform of the sequence x.
% The FRFT is computed using O(nlogn) operations.
%
% x       The sequence whose FRFT should be computed. Can be of odd or even
%         length. Must be a 1-D row vector.
% alpha	  The parameter alpha of the fractional Fourier transform
%
% Returns the aliased FRFT with parameter alpha of the sequence x.
% The fractional Fourier transform w of the sequence x (with parameter alpha) is defined by
%                   n/2-1
%       w(k) =       sum  x(u)*exp(-2*pi*i*k*u*alpha/N),  -n/2 <= k <= n/2-1, N=length(x).
%                   u=-n/2
%
% 
% Yoel Shkolnisky 18/12/02


function w=cfrftV2(x,alpha)
m=length(x);
%disp (strcat('FRFT LEN=',int2str(m)));
j=lowIdx(m):hiIdx(m);
j2= lowIdx(2*m):hiIdx(2*m);
E=i*pi*alpha;

paddedsize = pow2(nextpow2(3*m));
leftpad = fix((paddedsize-m+1)/2);
rightpad = paddedsize-m-leftpad; 
y=x.*exp(-E*j.^2/m);
y=[zeros(1,leftpad),y,zeros(1,rightpad)];

% compute the fractional Fourier transform not by padding to 3*m (as in the paper) but by
% padding to the next power of 2. This uses more memory, but since all FFTs use dyadic length
% it is much faster

z=zeros(1,paddedsize);
l=toUnaliasedIdx(-m,paddedsize);
z(l:l+length(j2)-1)=exp(E*j2.^2/m);

Y=cfft(y);
Z=cfft(z);
W=Y.*Z;
w=icfft(W);
w=w(toUnaliasedIdx(lowIdx(m),paddedsize):toUnaliasedIdx(hiIdx(m),paddedsize));
w=w.*exp(-E*j.^2/m);