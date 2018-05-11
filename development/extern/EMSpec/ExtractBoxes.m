function [images,pts0]=ExtractBoxes(imgname, boxname, n, imscale, displayOn,tMat)
% % function [images,pts]=ExtractBoxes(imgname, boxname, n, display,tMat)
% % Given a boxer particle coordinate file, create the images stack of n x n
% % images.  If display=1, the box locations are shown superimposed on the
% % image.
rotNum=0;
if nargin<6
    tMat=eye(3);
end;
itMat=inv(tMat);
if nargin<5
    displayOn=1;
end;
if nargin<4
    imscale=1;
end;
if nargin<3
    n=64;
end;

if nargin<2  % no filenames
    displayOn=1;
    [imgname pa]=uigetfile('*','Select an image file');
    cd(pa)
    [boxname pa]=uigetfile('*','Select a box coordinate file');
end;

if isnumeric(imgname)
    img=imgname;
else
    disp(['Reading image    ' imgname]);
    [img, pixA]=ReadEMFile(imgname);
    img=RemoveOutliers(img);  % *********
end;

disp(['Reading box file ' boxname]);
bfile=fopen(boxname);
if bfile<1
    error(['invalid file name ' boxname]);
end;
bc=textscan(bfile,'%n %n %n %n %*[^\n]');
fclose(bfile);
%%

n0=size(img);

img=rot90(single(img-mean(img(:))),rotNum); %!!!!!!!!
% Compute the centers of the boxes
% pts=[nx-bc{2}-bc{4}/2+1 bc{1}+bc{3}/2+1];
pts0=(imscale*([bc{1} bc{2}]+[bc{3} bc{4}]/2)+1)';  % This is what EMAN boxer does.
% pts=[nx-bc{2} bc{1}];
% pts=[bc{1} nx-bc{2}];
% pts=[bc{2} nx-bc{1}];
[nc,nim]=size(pts0);
limits=[min(pts0,[],2) max(pts0,[],2)]

disp(['Extracting ' num2str(nim) ' images']);
images=single(zeros(n,n,nim));

pts=AffineTransformPoints(pts0,n0,tMat);

for i=1:nim
images(:,:,i)=ExtractImage(img,round(pts(:,i)),n);
end;

if displayOn
    %
    % % Display the image and boxes
    w=3;  % box line width
mxi=max(img(:));
mni=min(img(:));
imgf=(img-mni)/(mxi-mni)*255;
%     imgf=single(imscale(img));

    [nb n2]=size(pts);

    box=zeros(n+2,n+2);

    box(1:w,:)=1; box(:,1:w)=1;
    box(n-w:n,:)=1; box(:,n-w:n)=1;
    m=255;
    boxi=1-box;
    imgd=Mask(imgf,pts,boxi,m*box);
    imacs(imgd);
end;

