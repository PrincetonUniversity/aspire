function [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
    cryo_parse_Relion_CTF_struct(CTFdata)

voltage=CTFdata.rlnVoltage;
DefocusU=CTFdata.rlnDefocusU/10; % Relion uses Angstrom. Convert to nm.

if isfield(CTFdata,'rlnDefocusV')
    DefocusV=CTFdata.rlnDefocusV/10; % Relion uses Angstrom. Convert to nm.
else
    DefocusV=DefocusU;
end

if isfield(CTFdata,'rlnDefocusAngle')
    DefocusAngle=CTFdata.rlnDefocusAngle*pi/180; % Convert to radians.
else
    DefocusAngle=0;
end

Cs=CTFdata.rlnSphericalAberration; % In mm, No conversion is needed.

if isfield(CTFdata,'rlnDetectorPixelSize')
    PS=CTFdata.rlnDetectorPixelSize; % In microns. Convert to Angstroms below.
    mag=CTFdata.rlnMagnification;
    pixA=PS*10^4/mag; % Convert pixel size on the detector in microns to spatial resoution in Angstroms.
elseif isfield(CTFdata,'pixA')
    pixA=CTFdata.pixA;
else
    error('Cannot get pixel size from CTF data');
end
A=CTFdata.rlnAmplitudeContrast;
