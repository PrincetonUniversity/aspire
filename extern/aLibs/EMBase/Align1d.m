function [v offs]=Align1d(in,ref)
% function [v offs]=Align1d(in,ref)
% Find the maximum cross-correlation between in and ref, and circularly
% shift in by the integer offs to best match the reference to yield the
% output vector v.

ccf=real(ifft(fft(in).*conj(fft(ref))));
[mxv i]=max(ccf);
offs=1-i;
v=circshift(in,offs);
