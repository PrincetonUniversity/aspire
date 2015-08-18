function mc=cryo_crop(m,n,isstack,fillval)
% CRYO_CROP     Crop proejctions
%
% mc=Crop(m,n,isstack,fillvalue)
%   Reduce the size of the 1d array, square or cube m by cropping (or
%   increase the size by padding with fillval, by default zero) to a final
%   size of n x n or n x n x n.  This is the analogue of Downsample, but
%   doesn't change magnification. 
%   If m is 2-dimensional and n is a vector, m is cropped to n=[nx ny].
%   The function handles odd and even-sized arrays correctly  The center of
%   an odd array is taken to be at (n+1)/2, and an even array is n/2+1.
%   If the flag isstack = 1 then a 3D array m is treated as a stack of 2D
%   images, and each image is cropped to n x n.
%   For 2D images, the input image doesn't have to be square.
%   The result is double if fillval is double; by default the result is
%   single.
%
% Written by Fred Sigworth. 
% Revised by Yoel Shkolnisky, August 2015.

if nargin<3
    isstack=0;
end;
if nargin<4
    fillval=single(0);  % Force a single output when padding.
end;

sz=size(m);
ndi=ndims(m);
if ndi==2 && any(sz==1)
    ndi=1;
end;

switch ndi
    case 1
        n1=numel(m);
        m=reshape(m,n1,1); % force a column vector
        ns=floor(n1/2)-floor(n/2);  % Shift term for scaling down.
        if ns>=0 % cropping down
            mc=m(ns+1:ns+n);
        else
            mc=fillval*ones(n,1);
            mc(1-ns:n1-ns)=m;
        end;
              
    case 2
        if numel(n)<2
            n(2)=n(1);
        end;
        nx=size(m,1);
        ny=size(m,2);
        nsx=floor(nx/2)-floor(n(1)/2);  % Shift term for scaling down.
        nsy=floor(ny/2)-floor(n(2)/2);
        if nsx>=0 % cropping down
            mc=m(nsx+1:nsx+n(1),nsy+1:nsy+n(2));
        else  % padding
            mc=fillval*ones(n);
            mc(1-nsx:nx-nsx,1-nsy:ny-nsy)=m;
        end;
        
    case 3 % m is 3D
        if isstack % a stack of 2D images
        if numel(n)<2
            n=n(1)*[1 1];
        end;
        ns=floor(sz(1:2)/2)-floor(n(1:2)/2);  % Shift term for scaling down.
            n(3)=sz(3);
            if all(ns>=0) % cropping down
                mc=m(ns(1)+1:ns(1)+n(1),ns(2)+1:ns(2)+n(2),:);
            elseif all(ns<=0)  % padding
                mc=fillval*ones(n);
                mc(1-ns(1):sz(1)-ns(1),1-ns(2):sz(2)-ns(2),:)=m;
            else
                error('Can''t crop and pad the same image');
            end;
        else  % not a stack
        if numel(n)<3
            n=n(1)*[1 1 1];
        end;
        ns=floor(sz/2)-floor(n/2);  % Shift term for scaling down.
            if all(ns>=0) % cropping down
                mc=m(ns(1)+1:ns(1)+n(1),ns(2)+1:ns(2)+n(2),ns(3)+1:ns(3)+n(3));
            elseif all(ns<=0)
                mc=fillval*ones(n);
                mc(1-ns(1):sz(1)-ns(1),1-ns(2):sz(2)-ns(2),1-ns(3):sz(3)-ns(3))=m;
            else
                error('Can''t crop and pad dimensions simultaneously');
            end;
        end;
    otherwise
        error(['Crop: dimension too large: ' num2str(ndi)]);
end;
