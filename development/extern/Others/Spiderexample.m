% read a SPIDER image and display it.
im = readSPIDERfile('img001.dat');

% the resulting matix is type double. To display using
% imshow or imtool, rescale to (0..1) with mat2gray :
imshow(mat2gray(im));
pause(2);
close

% display images from a stack using Matlab's montage:
s = readSPIDERfile('stk001.hcc');
t = mat2gray(s);

% For some reason, montage wants an extra singleton dimension.
% Create it by reshaping the array.
[ix iy iz] = size(t);
montage(reshape(t, [ix,iy,1,iz]));
pause(2);
close

% Extract some slices from a volume
v = readSPIDERfile('vol001.hcc');
[ix iy iz] = size(v);
k = round(ix/2);
imx =  squeeze(v(k,:,:));   % extract a slice from the x plane
imy =  squeeze(v(:,k,:));   % extract a slice from the y plane
imz =  squeeze(v(:,:,k));   % extract a slice from the z plane
% squeeze removes the extra singleton dimension, giving 2D images

% display them as subimages
subplot(1,3,1), subimage(mat2gray(imx))
subplot(1,3,2), subimage(mat2gray(imy))
subplot(1,3,3), subimage(mat2gray(imz))
pause(2);
close

% concatenate the images and write out as an image stack
C = cat(3,imx,imy,imz);
writeSPIDERfile('mystack.dat', C, 'stack')
disp('data written to mystack.dat')

% Display some slices through the volume
h = slice(double(v),k,k,k);
set(h,'FaceColor','interp','EdgeColor','none')
colormap gray
pause(3)
close

%----------------------------------------------
% document file functions

docfile = 'doc001.txt';
disp(sprintf('The contents of %s:', docfile))
type(docfile)
q = readSPIDERdoc(docfile) % read the entire file
x  = q(:,1);  % the first column
y1 = q(:,2);  % the second column
y2 = q(:,3);  % third column
plot(x,y1,x,y2)
pause(2)
close

% create some arrays and write them out
c1 = [1:10]       % c1, c2 are row vectors,use transpose (')
c2 = [1:10:100]   % to turn them into columns.
b2 = [c1' c2'];    % [ ] concatenates columns horizontally
outdoc = 'newdoc.txt';
writeSPIDERdoc(outdoc, b2, {'ones'; 'tens'});
disp(sprintf('data written to %s', outdoc))
type(outdoc)
