 function W = VesicleDensity(n,a,d,org)
% function W = VesicleDensity(n,a,d,org)
% Simple vesicle density calculation.
% Uses the same algorithm as VIDensity.m, using an integral over each pixel
% size.
% Compute the scattering density of a uniform spherical membrane of
% radius a, thickness d, center at org in an nxn image.
% The argument org is optional; org=[1 1] corresponds to the 
% lower-left corner of the image.
% Computed using the integral over pixels for greater accuracy.

% We subtract the density of a sphere a-d/2 in radius from a sphere of
% a+d/2 in radius.  The sphere density is computed in this way: 
% Each pixel is assigned a radius value r.  The value at that
% pixel is obtained as the integral of f(x)=sqrt(a^2-x^2) from r-1/2 to
% r+1/2.  In this way the resulting density W is smooth.

if nargin<4
    org=fctr(n);
end;

% % Simple calculation
% R=Radius(n,org);
% a1=a+d/2;
% a0=a-d/2;
% 
% W=2*real(sqrt(a1^2-R.^2)-sqrt(a0^2-R.^2));
% 
% return

% % fancy calculation
a1=a+d/2;
a0=a-d/2;
R=Radius(n,org);
Rp=R+0.5;
Rm=R-0.5;
W1=real(Rp.*sqrt(a1^2-Rp.^2)+a1^2*asin(Rp/a1)...
       -Rm.*sqrt(a1^2-Rm.^2)-a1^2*asin(Rm/a1));

W0=real(Rp.*sqrt(a0^2-Rp.^2)+a0^2*asin(Rp/a0)...
        -Rm.*sqrt(a0^2-Rm.^2)-a0^2*asin(Rm/a0));
W=(W1-W0);  % don't divide by 2 because we have 2 hemispheres.

% W=2*real(sqrt(a1^2-R.^2)-sqrt(a0^2-R.^2));
