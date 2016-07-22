function out=ExpandImage(in,nb)
% function out=ExpandImage(in,nb)
% Reverse of BinImage.  Presently works only with 2D image

if nb<=1
    out=in;
    return;
end;
nb=round(nb);

dims=ndims(in);
[nx ny nim]=size(in);
n=[nx ny];
dims=2;
if any(n==1)  % a 1d function
    dims=1;
    n=prod(n);
    in=reshape(in,n,nim);
end;
if numel(nb)<2  % if nb is a scalar, make it a vector
    nb=nb(1)*ones(1,2);
end;

if dims ~=2
    error('Only 2d images handled.');
end;

in=in';  % First, take the transpose
q=repmat(in(:)',nb(2),1);  % expand in y direction
qr=reshape(q,ny*nb(2),nx)'; % transpose back
out=repmat(qr(:)',nb(1),1); % expand in x direction
out=reshape(out,nx*nb(1),ny*nb(2));
