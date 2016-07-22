function h=ComplexCTF(n, pixA, Pars)
% function h=ComplexCTF(n, pixA, Pars)
% Compute the complex form of the contrast transfer function corresponding
% an n(1) x n(2) image with
% the sampling interval pixA in A/pixel.
%  Pars has the fields
%     lambda    (in A)
%     defocus   ( in um)
%     Cs        ( in mm)
%   and optionally
%     deltadef  (in um)
%     theta     (radians)
% The envelope constant B is ignored.
% alpha is also ignored.
% The result is returned in an nxn matrix with h(n/2+1,n/2+1) giving
% the zero-frequency amplitude.  Thus you must use ifftshift() on
% the result before performing a Fourier transform.
%
% The first zero (i.e. h is real) occurs at lambda*defocus*f0^2=1.
% e.g. when lambda=.025A, defocus=1um, then f0=1/16ÅA.
%
astig=0;
cuton=0;
lambda=Pars.lambda;
defocus=Pars.defocus;
Cs=Pars.Cs;
% alpha=Pars.alpha;
if isfield(Pars,'deltadef')  % we are handling astigmatism
    deltadef=Pars.deltadef;
    theta=Pars.theta;
    if deltadef ~= 0
        astig=1;
    end;
end;
if isfield(Pars,'cuton')
    cuton=P.cuton;
end;

if astig
    [r1 ang]=RadiusNorm(n,fctr(n));
    % we use the defocus formula from Henderson '86:
    df=defocus+deltadef*cos(2*(ang-theta));
else
    r1=RadiusNorm(n,fctr(n));
    df=defocus;
end;
r2=r1.^2;
f0 = 1./pixA;  % Spatial frequency unit (inverse Å)

k2=-df*pi*lambda*1e4*f0.^2;  % this may be a matrix
k4= pi/2*Cs*lambda^3*1e7*f0.^4;  % this is a scalar.
if Cs==0
    h=exp(1i*k2.*r2);
else
    h=exp(1i*(k2.*r2+k4*r2.*r2));
end

if cuton  % handle sharp cut-on of a phase plate.
    h=h.*(0.5+0.5*erf((r1*f0-cuton)*10/cuton));
end;
