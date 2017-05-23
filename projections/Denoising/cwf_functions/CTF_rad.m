function h = CTF_rad(n, pixA, lambda, defocus, Cs, B, alpha, r_val)
% f = CTF(n, pixA, lambda, defocus, Cs, B, alpha, r_val)% Cs term fixed.  fs 4 Apr 04 Struct option added fs 6 Jul 09
% Astig option added fs 18 Aug 09
% Tejal Bhamre: Oct 2015
% Evaluate CTF at radial points r_val
% No astigmatism, only radial function
r1=RadiusNorm(n,fctr(n));
df=defocus;
%Evaluate at sample points r_val only
r1=r_val;r2=r1.^2;
f0 = 1./pixA; 
 % Spatial frequency unit (inverse )
k2=-df*pi*lambda*1e4*f0.^2;  % this may be a matrix
k4= pi/2*Cs*lambda^3*1e7*f0.^4;  % this is a scalar.
kr=f0^2*B;  % B-factor
if Cs==0    
h=sin(k2.*r2-alpha).*exp(-kr*r2);
else    
h=sin(k2.*r2+k4*r2.*r2-alpha).*exp(-kr*r2);
end
