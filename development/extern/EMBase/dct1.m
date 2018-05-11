function y=dct1(x)
% 1D discrete cosine transform
% simple implementation

n0=numel(x);
x1=zeros(2*n0,1);
x1(1:n0)=x;  % make a padded copy
f1=fft(x1);
y=2*real(f1(1:n0));
y(1)=y(1)/2;
