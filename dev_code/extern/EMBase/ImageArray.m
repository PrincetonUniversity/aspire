function mb=ImageArray(m,norm, nd, border,nxi,nyi,bval)
% function mb=ImageArray(m, norm, nd, border, nxi,nyi,bval)
% Create a big array mb into which the stack of images m is tiled.
% All arguments except m are optional).
% nd is the number of points per image (after downsampling).
% Set norm=1 if you want each image to have its pixel variance normalized to 1.
% border is the border width between tiles;
% nxi and nyi are the number of tiles in the x and y directions.
% bval is the value in the border between tiles (default 3).

if ndims(m)>3 % already sorted into rows and columns
    [nx ny nxi nyi]=size(m);
    nim=nxi*nyi;
else
    [nx ny nim]=size(m);
end;

if nargin<2
    norm=0;
end;

if nargin<3
    nd=nx;
end;

if nargin<4
    border=1;
end;

nxx=max(1,ceil(sqrt(nim)));
if nargin<5 && ndims(m)<4
    nxi=nxx;
end;
if nargin<6 && ndims(m)<4
    nyi=ceil(nim/nxi);
end;

if nargin<7
    bval=1;
end;

mb=ones(nxi*nd+(nxi-1)*border,nyi*nd+(nyi-1)*border)*bval;
for iy=1:nyi
    y0=(iy-1)*(nd+border);
    for ix=1:nxi
        index=ix+nxi*(iy-1);
        if index <= nim
            x0=(ix-1)*(nd+border);
            mp=m(:,:,index);
            if nd<nx
                mp=Downsample(mp,nd);
            end;
            if norm
                mp=normalize(mp);
            end;
            for y=1:nd
                mb(x0+1:x0+nd,y0+y)=mp(:,y);
            end;
        end;
    end;
end;

end

function out=normalize(in)
me=mean(mean(in));
in=in-me;
var=sum(sum(in.^2))/numel(in);
out=in/sqrt(var);
end
