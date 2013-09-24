function P=FourierShift(n,sh)
% function P=FourierShift(n,sh)
% Compute the complex exponential for shifting an image.  n is the size of
% the image (scalar for a square image or a 2-element vector for a
% rectangular image), and can be odd.  If n is a 3-element vector, then a
% 3D array is produced. sh =[dx dy] or [dx dy dz] contains the shifts in
% pixels. Positive values are shifts up and to the right.  In the returned
% matrix, P(1,1) is zero frequency.
% If sh is nim x 2 in size then a stack of complex images is created.
%
% e.g. to shift by dx, dy pixels,
%   fm=fftn(m);
%   P=FourierShift(size(m),[dx dy]);
%   fsh=real(ifftn(fm.*P));
% or, in one line,
%   fsh=real(ifftn(fftn(m).*FourierShift(size(m,1),[dx dy])));
%
if numel(n)<2
    n(2)=n(1);
end;
if numel(sh)<4  % Too small to be a stack of shifts
    sh=sh(:)';  % force a row vector;
end;

ndims=numel(n);
nim=size(sh,1); % one row for each output image
switch ndims
    case 2
        P=zeros([n nim]);
        for i=1:nim
            nx=n(1);
            ny=n(2);
            [X,Y]=ndgrid((-floor(nx/2):floor((nx-1)/2))/nx,...
                (-floor(ny/2):floor((ny-1)/2))/ny);
            P(:,:,i)=exp((-1j*2*pi)*fftshift((sh(i,1)*X+sh(i,2)*Y)));
        end;
    case 3
        if nim>1
            error('Can''t handle a stack of 3D volumes');
        end;
        nx=n(1);
        ny=n(2);
        nz=n(3);
        [X,Y,Z]=ndgrid((-floor(nx/2):floor((nx-1)/2))/nx,...
            (-floor(ny/2):floor((ny-1)/2))/ny,...
            (-floor(nz/2):floor((nz-1)/2))/nz);
        P=exp((-1j*2*pi)*fftshift((sh(1)*X+sh(2)*Y+sh(3)*Z)));
    otherwise
        error('n must be a 1, 2 or 3-element vector');
end;
