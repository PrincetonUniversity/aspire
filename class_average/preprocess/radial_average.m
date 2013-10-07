function [ Pavg ,negmask] = radial_average( P )
%This function radially average the power spectrum
%      input:
%           P--original 2d power spectrum
%      output:
%           Pavg--1d power spectrum after radial averaging
%   Zhizhen Zhao June 2013

L=size(P, 1); %size of the power spectrum, accept odd dimension
N=floor(L/2);   
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
radial = zeros(ceil(N*sqrt(2))+1, 1);
for i=1:length(radial)
    radial(i)=mean(P(round(r)==i-1));
end;

Pavg = zeros(size(P));
for i=1:L
    for j = 1:L
        Pavg(i, j) = radial(round(r(i, j))+1);
    end;
end;

negmask=zeros(size(Pavg));
negmask(Pavg<0)=Pavg(Pavg<0);
Pavg(Pavg<0)=0;