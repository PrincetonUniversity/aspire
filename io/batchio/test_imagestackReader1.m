% Test the class imagestackReader
%
% The test script generates a stack of images where image i in the stack
% shows the number i. Then it loads the images from the stack, and the user
% should verify that the index of the images, which appears in the title of
% the figure matches the content of the image.
% See documentation in the script for more details instructions.
%
% Yoel Shkolnisky, May 2016.

% Test parameters
testfname='test.mrc'; % Name of the file containing test images.
sz=129;               % Size of each test image.  
N=35;               % Number of test image;

open_log(0);
% Create test data
log_message('Generating test data\n');
makeTestStack(testfname,sz,N);

% Create imagestack object
imstack=imagestackReader(testfname,5,'single',1);

% Read the first 5 images sequentially. No cache miss should happen (except
% for the first image. Make sure you get the correct images.
for k=1:5
    im=getImage(imstack,k);    
    imagesc(im);
    axis image; colormap(gray);
    title(sprintf('Verify this is image %d',k));
    log_message('Verify you see image %d',k);
    pause;
end


% Read images 5 to 10. No cache miss should happen (except
% for the first image. Make sure you get the correct images.
for k=5:10
    im=getImage(imstack,k);    
    imagesc(im);
    axis image; colormap(gray);
    title(sprintf('Verify this is image %d',k));
    log_message('Verify you see image %d',k);
    pause;
end


% Read images 25 and 35 (not in cache).
for k=[25 35]
    im=imstack.getImage(k);
    imagesc(im);
    axis image; colormap(gray);
    title(sprintf('Verify this is image %d',k));
    log_message('Verify you see image %d',k);
    pause;
end

% Read non-existing image
try
    im=imstack.getImage(N+1); % Should return an error.
catch
    log_message('Image %d returned an error since it does not exist. This is ok.',N+1);
end