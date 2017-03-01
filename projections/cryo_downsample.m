function [out, fx]=cryo_downsample(in,szout,stack,mask)
% CRYO_DOWNSAMPLE   down/up sample projections
%
% [out fx]=cryo_downsample(in,szout,stack,mask)
%   Use Fourier methods to change the sample interval and/or aspect ratio
%   of any dimensions of the input image 'in'. If the optional argument
%   stack=1 the last dimension of 'in' is interpreted as the index of each
%   image in the stack. The size argument szout is either a scalar or a
%   vector of the dimension of the output images.  Let the size of a stack
%   of 2D images 'in' be n1 x n1 x ni.  The size of the output (szout=n or
%   szout=[n n]) will be n x n x ni. The size argument szout can be chosen
%   to change the aspect ratio of the output; however the routine will not
%   allow one dimension to be scaled down and another scaled up.
%
%   If the optional mask argument is given, this is used as the
%   zero-centered Fourier mask for the resampling.  The size of mask should
%   be the same as the output image size. For example for downsampling an
%   n0 x n0 image with a 0.9 x nyquist filter, do the following:
%       msk=fuzzymask(n,2,.45*n,.05*n);
%       out=Downsample(in,n,0,msk);
%   The size of the mask must be the size of out. The optional fx output
%   argument is the padded or cropped, masked, FT of in, with zero
%   frequency at the origin. 
%
% Modified to operate on complex input: fs 15 Sep 2010.
% Modified to operate on rectangular 2D images fs Jan 2012
%
% Written by Fred Sigworth. 
% Revised by Yoel Shkolnisky, July 2015.


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
    mask=cryo_crop(mask,szout);  % scaling down: force it to the output size
end;
% elseif down
% else
%     mask=cryo_crop(mask,szin); % scaling up: use input size
% end;

% ns=(szin-szout)/2;  % shift
if ~copy
    if ~isa(in,'double');
    out=single(zeros([szout nim]));
    else
        out=zeros([szout nim]);
    end;
end;

showprogress=(nim>1000); % Determine if to show progress bar. Show progress 
                         % bar only for longer operations.

% Print progess bar only if processing more than 1000 images
if showprogress
    printProgressBarHeader;
end

switch ndim
    case 3
        if down  % scaling down
            for i=1:nim
                if showprogress
                    progressTicFor(i,nim);
                end
                x=fftshift(fftn(in(:,:,:,i)));
                fx=cryo_crop(x,szout).*mask;
                if ~copy
                    out(:,:,:,i)=ifftn(ifftshift(fx))*(prod(szout)/prod(szin));
                end;
            end;
        else      % scaling up
            for i=1:nim
                if showprogress
                    progressTicFor(i,nim);
                end
                x=fftshift(fftn(in(:,:,:,i)));               
                fx=cryo_crop(x,szout).*mask;
                out(:,:,:,i)=ifftn(ifftshift(fx))*(prod(szout)/prod(szin));
            end;
        end;
    case 2
        if down
            for i=1:nim
                if showprogress
                    progressTicFor(i,nim);
                end                
                fx=cryo_crop(fftshift(fftn(in(:,:,i))),szout).*mask;
                if ~copy
                    out(:,:,i)=ifftn(ifftshift(fx))*(prod(szout)./prod(szin));
                end;
            end;
        else   % scaling up
            for i=1:nim
                if showprogress
                    progressTicFor(i,nim);
                end                
                fx=cryo_crop(fftshift(fftn(in(:,:,i))),szout).*mask;
                out(:,:,i)=ifftn(ifftshift(fx))*(prod(szout)/prod(szin));
            end;
        end;
    case 1
        out=zeros([szout nim]);
        if down  % scaling down
            for i=1:nim
                if showprogress
                    progressTicFor(i,nim);
                end                
                fx=cryo_crop(fftshift(fft(in(:,i))),szout).*mask;
                if ~copy
                    out(:,i)=ifft(ifftshift(fx))*(prod(szout)/prod(szin));
                end;
            end;
        else   % scaling up
            for i=1:nim
                if showprogress
                    progressTicFor(i,nim);
                end                
                fx=cryo_crop(fftshift(fft(in(:,i))),szout).*mask;
                out(:,i)=ifft(ifftshift(fx))*(szout/szin);
            end;
        end;
        
end;
if isreal(in)
    out=real(out);
end;
fx=ifftshift(fx);  % shift back to the origin.