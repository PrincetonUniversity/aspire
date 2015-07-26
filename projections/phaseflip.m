function PFprojs=phaseflip(CTFdata,projs,prefix)
%
% Wrapper for backward compitability. 
%
% Yoel Shkolnisky, July 2015

warning('This function is depreacted. Call cryo_phaseflip instead.');
if nargin==2
    PFprojs=cryo_phaseflip(CTFdata,projs);
else
    PFprojs=cryo_phaseflip(CTFdata,projs,prefix);
end
