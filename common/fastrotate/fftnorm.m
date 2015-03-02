function z=fftnorm(signal);
%NORMALIZED DIRECT FFT
N=sqrt(size(signal,2));
signal=signal/N;
z=fft(signal);