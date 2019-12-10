function P=GetDefaultCTFPars(defocus,B,kV)
% function P=DefaultCTFPars(defocus,B,kV)
% Get default values for the structure P.  The following fields are set
% (parentheses indicate default values)
%   (kV=200)
%   P.lambda=EWavelength(kV)
%   P.defocus (1)
%   P.Cs=2
%   P.B (100+100*defocus)
%   P.alpha=.07

if nargin<1
    defocus=1;
end;
if nargin<2
    B=100+100*defocus;
end;
if nargin<3
    kV=200;
end;
P.lambda=EWavelength(kV);
P.defocus=defocus;
P.Cs=2;
P.B=B;
P.alpha=.07;

