function [aligned, pars, ccfs]=AlignTrans(in, ref, maxshift, oldpars)
% function [aligned, pars, ccfs]=AlignTrans(in, ref, maxshift, oldpars)
% Translational alignment of square images, making integer shifts.
% A stack of images 'in' is shifted so each aligns best with 'ref'; the
% output stack is 'aligned'.  maxshift=inf allows any size shift.
% The structure pars contains the alignment information and is updated
% from the (optional) old version.
% The ccfs output is the cross-correlation, 'fftshifted' to the center of
% the image.
% Note: at present the code doesn't support the pars calculation of shifts after
% rotations.
% Note: refs.shift is the shift of the input image needed to bring it into
% alignment with the reference.
% fs 27.6.04

[n n1 nim]=size(in);
cf=n/2+1;  % Center of F.T.

if nargin>2
    % Create a mask for the peak search.
    mask=zeros(n,n);
    lb=max(1,cf-maxshift);
    ub=min(n,cf+maxshift);
    mask(lb:ub,lb:ub)=1;
else
    mask=ones(n,n);
end;

fpre = conj(fft2( ref )); % We will use this to do correlation.

for i= 1:nim
    % Do the correlation
    fdat=fft2(in(:,:,i)); % FT of datum
    ycorr = ifftshift(real(ifft2(fdat.*fpre))); % cross-correlation of ypre and datum
    
    [maxval, is, js] = max2d(mask.*ycorr);
    pars.sx(i)=cf-is;  % zero means no shift
    pars.sy(i)=cf-js;  % zero means no shift
    if nargout>2
        ccfs(:,:,i)=ycorr;
    end;
end;

if nargin>3  % This needs to be fixed sometime to handle rot followed by trans.
    pars.sx=pars.sx+oldpars.sx;
    pars.sy=pars.sy+oldpars.sy;
    pars.rot=oldpars.rot;
else
    pars.rot=zeros(nim,1);
end;

% if nargout>1	% make output arrays
    for i=1:nim
        aligned(:,:,i) = shift( in(:,:,i), pars.sx(i), pars.sy(i) );
    end;
% end;
if nargout>2
    for i=1:nim
        ccfs(:,:,i)=ycorr;
    end;
end;