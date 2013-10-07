function [out finalmag nix]=DownsampleGeneral(in,nout,mag,maxOversample)
% function [out finalmag nix]=DownsampleGeneral(in,nout,mag,maxOversample)
% Perform 2D or 3D resampling but changing magnification by an arbitrary
% ratio.  Supports odd and even-sized square or cube inputs.  The ratio of
% approximated by the ratio noutx/ninx where ninx and noutx is a padded
% copies of the input and output, respectively.  maxOversample, which is
% by default 1.5 for large (>128 pixel) inputs and outputs, limits the
% maximum padding.  The finalmag value, the actual magnification ratio 
% which approximates mag, is returned.  Error are typically < 0.1%

% tol=9e-3;
ni=size(in,1);
ndim=sum(size(in)>1);

if nargin<3
    mag=nout/ni;
end;
if nargin<4
    nx=max(ni,nout/mag);
    if nx>128
        maxOversample=1.3;
    else
        maxOversample=256/(nx+16);
    end;
end;
% 

nos=nout:maxOversample*nout;
nis=round(nos/mag);
errs=(nos/mag-nis)./nis;  % relative errors
[val ind]=min(abs(errs));
% ind=find(abs(errs)<tol,1,'first');
% if numel(ind)<1
%     error('tolerance too small');
% end;

% subplot(212); plot(errs,'.-'); title(nis(ind));
% nix=nis(ind)
% nox=nos(ind)
finalmag=nos(ind)/nis(ind);

inp=Crop(in,nis(ind));  % pad the input
fout=Crop(fftshift(fftn(ifftshift(inp))),nos(ind));
out=Crop(fftshift(real(ifftn(ifftshift(fout)))),nout)*finalmag^ndim;

nix=nis(ind);