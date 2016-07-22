function [fts positions sds]=ComputeTiledFTs(mc,nu,nv,window)
% function [fts pos sds]=ComputeTiledFTs(mc,nu,nv,window)
% Given a large image mc, compute fts of tiles of size nu x nu which
% overlap by a minimum of nv pixels.  Before computing the fft, the
% windowed region is multiplied by window, if given; otherwise a circular
% window is used.  The mean of the windowed region is forced to zero, as is
% the edge of the window.  Return the omplex FTs, the positions of the 
% tiles, and the standard deviations of the tiles.  The vector
% positions(:,ix,iy)+1 gives the lower left pixel of the tile.


if nargin<3
    nv=0;  % default is no overlap.
end;
if nargin<4  % make the default window
    window=fuzzymask(nu,2,0.45*nu,0.1*nu);
end;
winsum=sum(window(:));

% Find out how many tiles we will use.
[nx ny]=size(mc);
tx=TileCoords(nx,nu,nv);  % no. of tiles in X direction
ty=TileCoords(ny,nu,nv);  % no. of tiles in Y direction

fts=complex(zeros(nu,nu,tx,ty));
positions=zeros(2,tx,ty);

if nargout>2
    sds=zeros(tx,ty);
end;

% loop over each tile.
for ix=1:tx
    x0=TileCoords(nx,nu,nv,ix);
    for iy=1:ty
        y0=TileCoords(ny,nu,nv,iy);
        positions(:,ix,iy)=[x0;y0];
        tm=mc(x0+1:x0+nu,y0+1:y0+nu);
        tm=tm.*window;
        tm=tm-sum(tm(:))*window/winsum;  % make the mean zero
        fts(:,:,ix,iy)=fftn(tm);
        if nargout>2
            sds(ix,iy)=sqrt(tm(:)'*tm(:))/nu;
        end;
    end;
end;
