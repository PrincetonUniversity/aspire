function [ proj ] = Phase_Flip( proj, defocus_group, c )
%This function is for phase flipping the images with give c
%   Detailed explanation goes here

% addpath /u/zhizhenz/matlab_cryo/common Need function cfft2.m and icfft2.m
L=size(proj, 1);
n=size(proj, 3);
K=size(c, 1);
l=floor(L/2);
k=ceil(K/2);
c=c(k-l:k+l, k-l:k+l, :);
for i=1:n
    proj(:, :, i)=icfft2(cfft2(proj(:, :, i)).*sign(c(:, :, defocus_group(i))));
end;

end

