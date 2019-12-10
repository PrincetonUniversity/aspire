function x=sect(m,dim)
% function x=sect(m,dim);
% Compute a central section of the image or volume m along the dimension dim.
% The output is then a central line or plane.
% if dim==0, a 3D input m is treated as a stack of nim images, and x is
% returned as an n x nim matrix of central lines along the first dimension.
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
        x=single(zeros(sz(1),nim));
        for i=1:nim
            x(:,i)=m(:,nct(2),i);
        end;
    case 3
        x=m(:,:,nct(3));
end;