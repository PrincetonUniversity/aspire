function [ Coeff ] = WF_FBSPCA( data, Mean, r_max, shifts, U, Freqs, W, filter_flag)
%This function computes the expansion coefficents
%   Input: data: images LxLxP
%          Mean: mean image
%          r_max: maximum radius
%          shifts: if region of interest is centered, put [0, 0]
%          U: principal components 
%          Freqs: The corresponding angular frequencies.
%          W: linear filter weights
%          filter_flag: 0---no Wiener filter
%                       1---Wiener Filter
%   Output: Coeff: expansion coefficents
%Zhizhen Zhao 09/01/2012

if nargin < 6
    filter_flag = 0;
end;
L=size(data, 1);
N=floor(L/2);
[x, y] = meshgrid(-N:N, -N:N);
r=sqrt((x-shifts(1)).^2+(y-shifts(2)).^2);
P=size(data, 3);

for i=1:P
    data(:, :, i)=data(:, :, i)-Mean;
end;

data=reshape(data, L^2, P);
data=data(r<=r_max, :);

Coeff=[U, conj(U(:, Freqs~=0))]\data;
Coeff=Coeff(1:length(Freqs), :);

if filter_flag == 1
    Coeff = diag(W)*Coeff;
end;
    
end

