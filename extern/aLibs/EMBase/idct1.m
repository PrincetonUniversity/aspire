function x=idct1(y)
% 1D inverse discrete cosine transform
% simple implementation

n0=numel(y);
y1=zeros(2*n0,1);
y1(1:n0)=y;  % make a padded copy
f1=ifft(y1);
x=2*real(f1(1:n0));
x(1)=x(1)/2;
