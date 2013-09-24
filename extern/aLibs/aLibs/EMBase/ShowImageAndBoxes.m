function mColor=ShowImageAndBoxes(m,ctrs,boxSize,boxLine,color,thresh)
% Display a grayscale image overlaid with boxes having the interior
% dimension boxSize and linewidth boxLine.  Color is a vector e.g.
% yellow = [1 1 0].
% ctrs is an nboxes x 2 array giving x,y pixel positions.
if nargin<6
    thresh=1e-3;
end;
if nargin<5
    color=[1 1 0];  % yellow
end;
n=size(m);
if nargin<4
    boxLine=round(n(1)/512);
end;
mScaled=single(imscale(m,1,thresh));  % scale to 0..1
innerBox=zeros(boxSize);
xBoxSize=boxSize+2*boxLine;
box=Crop(innerBox,xBoxSize,0,1);  % fill around it with ones
blank=zeros(xBoxSize);
blankMask=1-Mask(zeros(n),ctrs',blank,box);
drawMask=Mask(zeros(n),ctrs',blank,box);  % opaque boxes
% drawMask=Mask(zeros(n),ctrs,white,box);  % overlapping boxes
mColor=single(zeros([n 3]));
for i=1:3
    mColor(:,:,i)=mScaled.*blankMask+color(i)*drawMask;
end;
if nargout<1
    for i=1:3
        mColorR(:,:,i)=rot90(mColor(:,:,i));
    end;
    image(mColorR);
end;
