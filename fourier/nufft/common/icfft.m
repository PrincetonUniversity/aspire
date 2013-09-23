% function y=icfft(x)
%
% Aliased inverse Fourier transform (IFFT) of the sequence x.
% The FFT is computed using O(nlogn) operations.
%
% x    The sequence whose IFFT should be computed. 
%      Can be of odd or even length. Must be a 1D vector.
%
% Returns the aliased IFFT of the sequence x.
% 
% Yoel Shkolnisky 22/10/01

function y=icfft(x)
%y = fftshift1d(ifft(ifftshift1d(x)));

m=length(x);

if (mod(m,2)==1)
    pf=floor(m/2);
    pc=ceil(m/2);
else
    pf=m/2;
    pc=pf;
end

ix=x([pf+1:m 1:pf]);
ify=ifft(ix);
y=ify([pc+1:m 1:pc]);

% Revision Record
% 15/1/03	Yoel Shkolnisky		Use fftshift1d instead of fftshift
% 19/03/07  Yoel Shkolnisky     Remove calls to fftshift1d and ifftshift1d
