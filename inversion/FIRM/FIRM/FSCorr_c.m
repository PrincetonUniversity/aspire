function [c,f1,f2]=FSCorr_c(m1,m2,degree)
% c=FSCorr(m1,m2)
% Compute the fourier shell correlation between the 3D maps m1 and m2 in the missing cone,
% which must be n x n x n in size with the same n, assumed to be even.
%
% The FSC is defined as
%           \sum{F1 .* conj(F2)}
% c(i) = -----------------------------
%        sqrt(\sum|F1|^2 * \sum|F2|^2)
%
% Where F1 and F2 are the Fourier components at the given spatial frequency i.
% i ranges from 1 to n/2-1, times the unit frequency 1/(n*res) where res
% is the pixel size.

% First, construct the radius values for defining the shells.
[n n1 n2]=size(m1);
[x y z] = meshgrid(-n/2:n/2-1);  % zero at element n/2+1.
R = sqrt(x.^2 + y.^2 + z.^2);
eps = 1e-4;
% Creat a mask according to the missing cone
mask=zeros(n,n1,n2);
mask(abs(z)./sqrt(x.^2+y.^2)>tand(degree))=1;

% Fourier-transform the maps
f1=fftshift(fftn(ifftshift(m1)));
f2=fftshift(fftn(ifftshift(m2)));

f1=f1.*mask;
f2=f2.*mask;

% Perform the sums.
d0=R<0.5 + eps;  % Sphere at origin
c=zeros(n/2-1,1);
for i=1:n/2-1
	d1=R<0.5+i+eps;
	ring=d1-d0;
% 	nr=sum(sum(sum(ring)))
	r1=ring .* f1;
	r2=ring .* f2;
	num=real(sum(sum(sum(r1.*conj(r2)))));
	den=sqrt(sum(sum(sum(abs(r1).^2)))*sum(sum(sum(abs(r2).^2))));
	c(i)=num/den;
	d0=d1;
end;