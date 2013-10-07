function out=shiftf(in,delta)
% function out=shiftf(in,delta)
% Analogous to circshift, but uses the fft to  
% perform fractional shifts of the image in.
% delta = [dx dy] is the shift.
[nx ny]=size(in);
[X,Y]=ndgrid((-nx/2:nx/2-1)/nx,(-ny/2:ny/2-1)/ny);
P=exp(1j*2*pi*fftshift((delta(1)*X+delta(2)*Y)));
out=real(ifftn(fftn(in).*P));
end
