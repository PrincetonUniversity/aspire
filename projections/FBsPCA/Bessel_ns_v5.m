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

ns = unique(n);

[r_unique, ~, r_I] = unique(r(:));

for k = 1:numel(ns)
	nk = ns(k);
	Y = exp(i*nk*theta(:));

	mask = find(n==nk);

	r0 = r_unique*R_ns(mask(:))';
	F = besselj(nk, r0);
	F = F(r_I,:);

	phik = bsxfun(@times, Y, F);
	phik = bsxfun(@times, phik, 1./(N*sqrt(pi)*abs(besselj(nk+1, R_ns(mask(:))'))));

	phik = reshape(phik, [size(theta) numel(mask)]);

	Phi_ns(:,:,mask) = phik;
end

Phi_ns=reshape(Phi_ns, (2*N+1)^2, size(B, 1));
Phi_ns=Phi_ns(r<=1, :);

end
 
