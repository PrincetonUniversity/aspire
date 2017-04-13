function [vol,outstackname]=cryo_orient_projections_auxpreprocess_outofcore(vol,instackname)
%
% See cryo_orient_projections_auxpreprocess.m
%
% Yoel Shkolnisky, August 2015.

szvol=size(vol);

% Resample reference volume to the size of the projections
instack=imagestackReader(instackname);
szprojs=instack.dim; % Check that projections are square.
if szprojs(1)~=szprojs(2)
    error('Projections to orient must be square.');
end

if any(szvol-szvol(1)) % Check that the volume has all dimensions equal.
    error('Volume must have all dimensions equal');
end

% Mask and filter the volume
n=size(vol,1);
vol=vol.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));
%vol=GaussFilt(vol,0.3);

% Mask and filter the projections
tempstack1=tempmrcname;
outstack=imagestackWriter(tempstack1,instack.dim(3));
for k=1:instack.dim(3)
    p=instack.getImage(k);
    p=p.*fuzzymask(n,2,floor(0.45*n),floor(0.05*n));
    %p=GaussFilt(p,0.3);
    outstack.append(p);
end
outstack.close;

tempstack2=tempmrcname;
cryo_globalphaseflip_outofcore(tempstack1,tempstack2);
outstackname=tempstack2;