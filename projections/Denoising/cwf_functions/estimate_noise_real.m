function [ noise_variance ] = estimate_noise_real( data )
%% Estimating noise in real space images to give as input to ASPIRE, since it needs to be estimated before frequency thresholding
% but the input to ASPIRE is already frequency thresholded.
% Tejal Oct 2015

P=size(data, 3);
L=size(data, 1);
N=floor(L/2)-1;
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
r_max=N;
%finish set up paramters

%Estimate noise variance from the corners
test=reshape(data, L^2, P);
test=test(r>r_max, :);
noise_variance=var(test(:));


end

