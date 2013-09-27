function P=gridMakeNullFT(n, dims, mode)
% function P=gridMakeNullFT(n, dims, mode)
% Create the data structure for the gridding routines, assuming a real
% input structure of size n and dims from 1 to 3.
% The returned PadFT array is zeros.
%
% P is a data structure containing
%   n the original size
%   padfactor
%   np the initial padded size (np=padfactor*n)
%   np1 the size of the padded ft (larger than np)
%   sp the origin shift between n and np
%   sp1 the origin shift between np and np1
%   np1ctr the center of the ft after final padding
%   PadFT the final padded FT (each dimension np1), in CAS format.

if nargin<3
    mode='grid';
end;

padfactor=gridPadFactor(mode);

np=padfactor*n;  % Padding to do to original real image
np1=np+16; % Extra padding for ft of image
sp=ceil(np/2)-n/2; % shift of origin
sp1=ceil(np1/2)-np/2;  % shift of origin going from np to np1
np1ctr=n/2+1+sp+sp1;

x0=sp+1;
x1=sp+n;
y0=sp1+1;
y1=sp1+np;

% simply make zero arrays.
switch dims
    case 1
        P.PadFT=zeros(np1,1);
    case 2
        P.PadFT=zeros(np1,np1);
    case 3
        P.PadFT=zeros(np1,np1,np1);
    otherwise
        error(['Invalid dims argument: ',num2str(dims)]);
end;

P.n=n;
P.padfactor=padfactor;
P.np=np;
P.np1=np1;
P.sp=sp;
P.sp1=sp1;
P.np1ctr=np1ctr;
