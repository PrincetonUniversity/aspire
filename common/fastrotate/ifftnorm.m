function z=ifftnorm(signal);
%NORMALIZED INVERSE FFT
N=sqrt(size(signal,2));
signal=signal*N;
z=ifft(signal);