function [x0 x1 u1]=TileCoords(nx,nu,nv,index)
% function [x0 x1 u1]=TileCoords(nx,nu,nv,index)
% Get coordinates for tiling a large image along one dimension.
% There are two ways to call the function.
% nt=TileCoords(nx,nu,nv)  % Determine the number of tiles
% [x0 x1 u1]=TileCoords(nx,nu,nv,index)  % Compute the starting coordinate
% for a tile in a large input domain. If the function is called with
% index=0, it returns nt, the number of tiles. The input domain coordinate
% is denoted by x, total size is nx; the tile coordinate is u, with total
% size nu.  The minimum border (total overlap) between adjacent tiles is
% nv.
% Given an index [1..nt], the the returned values are the start location in
% the input domain x0 (x0+1 corresponds to u=1 in the tile), the location x1
% where data from the tile will be kept (x1 is approximately x0+nv/2), and
% the tile location u1 that corresponds to x1.

% Sample code:
% 
% m0=imread('dm2dwc01242c.tif','tif');
% [nx ny]=size(m0);
% m1=m0*0;
% nu=1024;
% nv=256;
% tx=TileCoords(nx,nu,nv);  % no. of tiles in X direction
% ty=TileCoords(ny,nu,nv);  % no. of tiles in Y direction
% 
% for ix=1:tx
%     [x0 x1 u1]=TileCoords(nx,nu,nv,ix);
% 
%     for iy=1:ty
%         [y0 y1 v1]=TileCoords(ny,nu,nv,iy);
% 
%         tm=double(m0(x0+1:x0+nu,y0+1:y0+nu))+offs;
%         tfilt=GaussFilt(tm,0.1);
%         m1(x1:x0+nu,y1:y0+nu)=uint8(tfilt(u1:nu,v1:nu));
%     end;
% end;

if nargin<4
    index=0;
end;

nt=max(1,ceil((nx-nv)/(nu-nv)));
if index==0
    x0=nt;  % Return nt as the first output argument.
    return
end;

if nx<=nu  % Trivial case: no tiling to do.
    x0=1;
    x1=1;
    u1=1;
    return;
end;

nv1=(nt*nu-nx)/(nt-1);  % Actual nv to use.
if index==1
    nv1=0;
end;

x0=(index-1)*(nu-nv1);  % starting copy location
x0=min(x0,nx-nu);  % Make sure we don't go off the end.
x0=ceil(x0);
x1=x0+1+round(nv1/2);
u1=x1-x0;

end
