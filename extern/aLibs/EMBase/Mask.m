function out=Mask(in,cpoint,mul,add)
% function out=Mask(in,cpoint,mul,add);
% Modify the vicinity of the point 'cpoint' of the input matrix, or the n
% points described by a 2 x n matrix.
% fs 19 Mar 03; changed for multiple points Feb 08.
% The matrix mul is centered on cpoint=[i j] and multiplied point-by-point with
% the corresponding elements of in.  Then the (optional) matrix add is
% added.  mul and add are assumed to be the same size.  If their dimensions
% n are odd, then their center (index (n+1)/2) is placed on (i,j); if even,
% the index n/2+1 is placed on (i,j), and there will be one more point below
% and to the left of the center than above and to the right (in cart.
% coordinates).
% Clipping is performed so that cpoint may be on the edge (or outside the
% bounds) of in.
% Related function: see ExtractImage.m

if nargin < 4
    add=0;
end;

szm=size(mul)';
lower=floor(szm/2); upper=szm-lower-1;
% For example, mapping 3x3 into (10,10):
% lower=1, upper = 1, use points i-1:i+1, center on point 2
% e.g. mapping 4x4
% lower=2, upper =1, use points i-2:i+1, center on point 3.
n=numel(cpoint);
cpoint=round(cpoint);
npoints=1;
cpts=cpoint;

out=in;

if n>2
    [n npoints]=size(cpoint);
    cpts=cpoint;
    cpoint=zeros(2,1);
elseif n==2
    cpts=cpoint(:);
elseif n<2
    return
end;


for i=1:npoints
    cpoint=cpts(:,i);

    % Clip the coordinates.
    inl=max([1 1]',cpoint-lower);  % lower left corner on the input
    inu=min(size(in)',cpoint+upper);
    rgl=inl-cpoint+lower+1; % lower left corner of the mask
    rgu=inu-cpoint-upper+szm;

    region=zeros(szm');
    region(rgl(1):rgu(1),rgl(2):rgu(2))=out(inl(1):inu(1),inl(2):inu(2));

    region=region.*mul+add;

    out(inl(1):inu(1),inl(2):inu(2))=region(rgl(1):rgu(1),rgl(2):rgu(2));
end;
