function [proj]=low_pass_filter(data)

P=size(data, 3);
L=size(data, 1);
k=floor(L/3);
proj=zeros(size(data));
for i=1:P
    pf=cfft2(data(:, :, i));
    pf2=zeros(L);
    pf2(k+1:L-k, k+1:L-k)=pf(k+1:L-k, k+1:L-k);
    proj(:, :, i)=icfft2(pf2);
end;