function x=sectr(m,dim)
% function x=sectr(m,dim)
% returns a 1D section of the 2D or 3D array m.
% This is like the sect() function except that the section runs only
% from the center to the maximum index along dimension dim.
% if dim==0, a 3D input m is treated as a stack of nim images, and x is
% returned as an n/2 x nim matrix of central lines along the first dimension.
% if nargin<2
%     dim=1;
% end;
% m=squeeze(shiftdim(m,dim-1));
% sz=size(m);
% nd=ndims(m);
% nct=ceil(sz/2+1);  % center of each dimension
% switch nd
%     case 2
%         x=m(nct(1):sz(1),nct(2));
%     case 3
%         x=m(nct(1):sz(1),nct(2),nct(3));
% end;
% 
% 
if nargin<2
    dim=0;
end;
if dim>1
    dim=dim-1;
    m=squeeze(shiftdim(m,dim));
end;
sz=size(m);
nd=ndims(m);
if dim<1 && nd>2
    nim=sz(nd);
    nd=nd-1;
    sz=sz(1:nd);
else
    nim=1;
end;
nct=ceil((sz+1)/2);  % center of each dimension
switch nd
    case 2
        x=zeros(ceil(sz(1)/2),nim);
        for i=1:nim
            x(:,i)=m(nct(1):sz(1),nct(2),i);
        end;
    case 3
        x=m(nct(1):sz(1),nct(2),nct(3));
end;