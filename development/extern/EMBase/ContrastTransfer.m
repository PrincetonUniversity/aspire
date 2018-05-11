function [c chi]=ContrastTransfer(s, lambda, defocus, Cs, B, alpha)
% function [c chi]=ContrastTransfer(s, lambda, defocus, Cs, B, alpha)
% function [c chi]=ContrastTransfer(s, CPars)
% Compute the contrast transfer function corresponding
% to the given spatial frequency s (in A^-1); lambda is in A, Cs is in mm,
% B in A^2 and alpha in radians.  Defocus is in microns.
% The returned chi variable is -1/pi times the argument to the sin function.
% Hence CTF=sin(-pi*chi).*exp(-B*s.^2). this form is convenient because
% chi=-1 at the first zero, -2 at the second, etc. as long as Cs effects
% are negligible (otherwise it's non-monotonic).
% In the alternate form, CPars is a structure containing those parameters.
% CPars.lambda
% CPars.defocus
% CPars.Cs
% CPars.B
% CPars.alpha
%   and, optionally, to handle astigmatism,
% CPars.deltadef
% CPars.theta
% The frequency s can be of any dimension.
% Note that astigmatism makes sense only for 2D frequency values.
% For 1D frequency arrays, the values along the x-axis (ang=0) are taken.

if isstruct(lambda)
    P=lambda;
    lambda=P.lambda;
    defocus=P.defocus;
    Cs=P.Cs;
    B=P.B;
    alpha=P.alpha;
    if isfield(P,'cuton')
        cuton=P.cuton;
    else
        cuton=0;
    end;
    if isfield(P,'deltadef')  % we are handling astigmatism
        deltadef=P.deltadef;
        theta=P.theta;
    else
        deltadef=0;
    end;
end;

if deltadef~=0 && ndims(s)<3  % handle astigmatism
    if all(size(s)>1)         % Truly 2-dimensional
        ang=atan2(s(2),s(1));
    else
        ang=0;
    end;
    defocus=defocus+deltadef*cos(2*(ang-theta));
end;
s2=s.^2;
chi=-1e4*lambda.*defocus.*s2+Cs.*lambda.^3.*5e6.*s2.^2-alpha/pi;

c=sin(pi*chi).*exp(-B*s2);
if cuton  % handle sharp cut-on of a phase plate.
    c=c.*(0.5+0.5*erf((abs(s)-cuton)*10/cuton));
end;
