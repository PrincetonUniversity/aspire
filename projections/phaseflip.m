function PFprojs=phaseflip(CTFdata,projs,prefix)
%
% Wrapper for backward compitability. 
%
% Yoel Shkolnisky, July 2015

warning('This function is depreacted. Call cryo_phaseflip instead.');
PFprojs=cryo_phaseflip(CTFdata,projs,prefix);
