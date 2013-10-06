function P=gridMakePaddedFT(m, mode, precomp)
% function P=gridMakePaddedFT(m, mode, precomp)
% Take a band-limited, floated image or volume m and compute its interpolated and
% padded Fourier transform for use with sinc or gridding interpolation.
% For sinc interpolation (no gridding), use padfactor = 2;
% for gridding use padfactor=1.25.
% The argument precomp is either 1 (for no pre-compensation) or a 1D array created by
% gridMakePrecomp to compensate for the roll-off of the interpolation
% kernel.  It is expanded to the proper dimension and multiplies m before
% tranforming.
% The output Fourier transform is presented as a real array of CAS values.
% Its dimension is determined from the dimension of m, either 2 or 3.
%
% P is a data structure containing
%   np the initial padded size (np=padfactor*n)
%   np1 the size of the padded ft (larger than np)
%   sp the origin shift between n and np
%   sp1 the origin shift between np and np1
%   np1ctr the center of the ft after final padding
%   PadFT the final padded FT of the input image m, in CAS format.

if nargin<3
    precomp=1;
end;
if nargin<2
    mode='grid';
end;

padfactor=gridPadFactor(mode);
sizes=size(m);
n=sizes(1);
dimension=ndims(m);
if sizes(2)==1
    dimension=1;
end;

np=padfactor*n;  % Padding to do to original real image
np1=np+16; % Extra padding for ft of image
sp=ceil(np/2)-n/2; % shift of origin
sp1=ceil(np1/2)-np/2;  % shift of origin going from np to np1
np1ctr=n/2+1+sp+sp1;

x0=sp+1;
x1=sp+n;
y0=sp1+1;
y1=sp1+np;

switch dimension
    case 1
        mp=zeros(np,1);
        mp(x0:x1)=m;
        mp=mp.*precomp;
        fmp=fftshift(fft(fftshift(mp)));
        fmp1=zeros(np1,1);
        fmp1(y0:y1)=fmp; % padded ft.
        fmp1=fmp1.*fuzzymask(np1,1,np/2+2,3,np1ctr); % mask outside frequencies.

    case 2
        % Create the padded realspace original
        mp=zeros(np,np);
        mp(x0:x1,x0:x1)=m;

% whos mp
% whos precomp
mp=mp.*kron(precomp,precomp');  % perform the compensation

        fmp=fftshift(fftn(fftshift(mp)));  % Get the ft, which is now fine-grained.

        % Pad the ft to avoid problems in the interpolation
        fmp1=zeros(np1,np1);
        fmp1(y0:y1,y0:y1)=fmp; % padded ft.
        fmp1=fmp1.*fuzzymask(np1,2,np/2+2,3,[np1ctr,np1ctr]); % mask outside frequencies.

    case 3
        % Make the expanded transform
        mp=zeros(np,np,np);
        mp(x0:x1,x0:x1,x0:x1)=m;

        % pre-compensate
        if numel(precomp)>1
            if numel(precomp)~=np
                error(['precomp is wrong size: ' num2str(numel(precomp)) ', np= ' num2str(np)]);
            end;
            mp=mp.*reshape(kron(precomp,kron(precomp,precomp)),np,np,np);  % perform the 3D compensation
        end;
        fmp=fftshift(fftn(fftshift(mp)));  % Get the ft, which is now fine-grained.
        % Pad the ft to avoid problems in the interpolation
        fmp1=zeros(np1,np1,np1);
        fmp1(y0:y1,y0:y1,y0:y1)=fmp; % padded ft.
        fmp1=fmp1.*fuzzymask(np1,3,np/2+2,3,[np1ctr,np1ctr,np1ctr]); % mask outside frequencies.
end;
P.PadFT=ToCAS(fmp1);

P.n=n;
P.np=np;
P.np1=np1;
P.sp=sp;
P.sp1=sp1;
P.np1ctr=np1ctr;
P.padfactor=padfactor;
