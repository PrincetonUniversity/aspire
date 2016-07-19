function [ freqs ] = pft_freqs(r, n_theta)
%%% Description
%The function 
%       1. Computes the frequencies <k, r> for polar Fourier transform
%
%Input: 
%       1. r: maximum angular frequency k_max
%       2. n_theta: number of samples on concentric circles
%Output:
%       1. freqs: list of frequencies <k, r>.
%Zhizhen Zhao 02/2015

n_r = length(r);
dtheta=2*pi/n_theta;

freqs=zeros(n_r*n_theta,2); % sampling points in the Fourier domain
for j=1:n_theta
      freqs((j-1)*n_r+1:j*n_r,:)=[r*sin((j-1)*dtheta),...
            r*cos((j-1)*dtheta)];
end
