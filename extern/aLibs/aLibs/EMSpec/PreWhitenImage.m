function out=PreWhitenImage(in)
% Given an unbinned image read from our Ultrascan camera at 200 keV, operate with an
% inverse filter to pre-whiten the noise.
n=size(in,1);
if n ~= 4096
    disp('PreWhitenImage warning: image size is not 4k');
end;
f=Radius(n)/n;  % normalized frequency
out=real(ifftn(fftn(in).*fftshift(CCDPreWhite(f))));
