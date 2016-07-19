% Test the class imagestackReader
%
% Test reading a bulk of images.
% The test script generates a stack of images where image i in the stack
% shows the number i. Then it loads the several images at once from the
% stack, and the user should verify that the indices of the images, which
% appear in the title of the figure match the content of the images.
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
imstack=imagestackReader(testfname,5);

idx=randperm(N);
idx=idx(1:16);
images=imstack.getImage(idx);    
viewstack(images,4,4);
title(sprintf('Verify that you see images \n %s',num2str(idx)));