function [out fx]=Downsample(in,szout,stack,mask)
% function [out fx]=Downsample(in,szout,stack,mask)
% Use Fourier methods to change the sample interval and/or aspect ratio of
% the dim=1,2 or 3 image in. If the optional argument stack=1 the (dim+1)
% dimension is interpreted to be the number of images in the stack. The
% size argument szout is either a scalar or a vector of the size of the
% output images.  Let the size of a stack of 2D images 'in' be n1 x n1 x
% ni.  The size of output (szout=n or szout=[n n]) will be n x n x ni.
%  The size argument szout can be chosen to change the aspect ratio of the
% output; however the routine will not allow one dimension to be scaled
% down and another scaled up.
%
% If the optional mask argument is given, this is used as the zero-centered
% Fourier mask for the resampling.  The size of mask should be the same as
% the output image size. For example for downsampling an n0 x n0 image with
% a 0.9 x nyquist filter, do the following:
%   msk=fuzzymask(n,2,.45*n,.05*n);
%   out=Downsample(in,n,0,msk);
% The size of the mask must be the size of out. The optional fx
% output argument is the padded or cropped, masked, FT of in, with zero frequency
% at the origin. 

% Modified to operate on complex input: fs 15 Sep 2010.
% Modified to operate on rectangular 2D images fs Jan 2012

if nargin<3
    stack=0;
end;

if isa(in,'integer')
    in=single(in);
    in=in-mean(in(:));
end;

nim=1;
ndim=sum(size(in)>1);  % number of non-singleton dimensions
if stack
    nim=size(in,ndim);
    ndim=ndim-1;
end;

szin=size(in);
szin=szin(1:ndim);

szout=szout(:)';  % force a row vector of size ndim
if numel(szout)<ndim
    szout=szout(1)*ones(1,ndim);
end;

if ndim==1
    in=in(:);  % force input to be a column vector in the 1d case
end;

copy=0;
if all(szout==szin)  % no change in size
    out=in;
    if nargout<2
        return
    else  % have to compute fx
        copy=1;
    end;
    
end;

down=all(szout<=szin);  % scaling down
% Make sure we aren't scaling down and up at the same time.
if ~down && any(szout<szin)  % Not all scaling up
    error('Incompatible dimension change');
end;

if nargin<4
    mask=1;
else
    mask=Crop(mask,szout);  % scaling down: force it to the output size
end;
% elseif down
% else
%     mask=Crop(mask,szin); % scaling up: use input size
% end;

% ns=(szin-szout)/2;  % shift
if ~copy
    if ~isa(in,'double');
    out=single(zeros([szout nim]));
    else
        out=zeros([szout nim]);
    end;
end;

switch ndim
    case 3
        if down  % scaling down
            for i=1:nim
                x=fftshift(fftn(in(:,:,:,i)));
                fx=Crop(x,szout).*mask;
                if ~copy
                    out(:,:,:,i)=ifftn(ifftshift(fx))*(prod(szout)/prod(szin));
                end;
            end;
        else      % scaling up
            for i=1:nim
                x=fftshift(fftn(in(:,:,:,i)));
                
                fx=Crop(x,szout).*mask;
                out(:,:,:,i)=ifftn(ifftshift(fx))*(prod(szout)/prod(szin));
            end;
        end;
    case 2
        if down
            for i=1:nim
                fx=Crop(fftshift(fftn(in(:,:,i))),szout).*mask;
                if ~copy
                    out(:,:,i)=ifftn(ifftshift(fx))*(prod(szout)./prod(szin));
                end;
            end;
        else   % scaling up
            for i=1:nim
                fx=Crop(fftshift(fftn(in(:,:,i))),szout).*mask;
                out(:,:,i)=ifftn(ifftshift(fx))*(prod(szout)/prod(szin));
            end;
        end;
    case 1
        out=zeros([szout nim]);
        if down  % scaling down
            for i=1:nim
                fx=Crop(fftshift(fft(in(:,i))),szout).*mask;
                if ~copy
                    out(:,i)=ifft(ifftshift(fx))*(prod(szout)/prod(szin));
                end;
            end;
        else   % scaling up
            for i=1:nim
                fx=Crop(fftshift(fft(in(:,i))),szout).*mask;
                out(:,i)=ifft(ifftshift(fx))*(szout/szin);
            end;
        end;
        
end;
if isreal(in)
    out=real(out);
end;
fx=ifftshift(fx);  % shift back to the origin.