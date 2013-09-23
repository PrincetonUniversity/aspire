function [ Phi_ns, n, s, R_ns ]=Bessel_ns_v5( N )
%This function generates the Fourier-Bessel Basis with positive angular frequencies
%up to 179
% Table Bessel180 gives the R_ns (3rd column) for n (1st column) and s (2nd column)
%   Input: N maximum radius
%   Output: Phi_ns: Fourier Bessel Basis with order n and s
%           n: angular frequencies
%           s: radial frequencies
%   Added Normalizing factor
% Zhizhen Zhao 2012 August

load Bessel180.mat %This is a table for Bessel Zeros
% if cutoff>2*N, cutoff = 2*N; end;

B = Bessel180;
clear Bessel180;

B = B(B(:, 3)<pi*N, :);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
theta=atan2(y, x);
r=r/N;
n = B(:, 1);
s = B(:, 2);
R_ns = B(:, 3);
Phi_ns=zeros(2*N+1, 2*N+1, size(B, 1));

%Go concurrent
ps=matlabpool('size');
if ps==0
    matlabpool open
end

parfor i=1:size(B, 1)
    r0=r*R_ns(i);
    [ F ]=besselj(n(i), r0); %this is the radial function
    Y=exp(sqrt(-1)*n(i)*theta);
    Phi_ns(:, :, i)=1/(N*sqrt(pi)*abs(besselj(n(i)+1, R_ns(i))))*F.*Y;
end;

Phi_ns=reshape(Phi_ns, (2*N+1)^2, size(B, 1));
Phi_ns=Phi_ns(r<=1, :);

end
 