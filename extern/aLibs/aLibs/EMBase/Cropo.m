function mc=Cropo(m,n,isstack,fillval)
% function mc=Cropo(m,n,isstack,fillval)
% Crop centered on origin.  This is the same as Crop,
% but assumes that m(1,1) is the 'center' of the image to be cropped.
% Equivalent to ifftshift(Crop(fftshift(m),n)), except it can handle stacks
% of 2d images too.
% For example, to downsample an image to nxn in size, you can do this:
% dm=real(ifftn(Cropo(fftn(m),n)));

if nargin<3
    isstack=0;
end;
if nargin<4
    fillval=single(0);  % Force a single output when padding.
end;

ndi=ndims(m);

if ndi==3 && isstack
   ni=size(m,3);
   ndi=2;
   mc=m;
   for i=1:ni
       mc(:,:,i)=ifftshift(Crop(fftshift(m(:,:,i)),n,0,fillval));
   end;
else
    mc=ifftshift(Crop(fftshift(m),n,0,fillval));
end;
