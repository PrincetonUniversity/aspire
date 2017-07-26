function [ noise_variance ] = estimate_noise_v6( data )
% Estimating noise from corners of images
% Tejal Oct 2015

P=size(data, 3);
L=size(data, 1);
N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
r_max=N;

%Estimate noise variance from the corners
test=reshape(data, L^2, P);
test=test(r>r_max, :);
noise_variance=var(test(:));


end

