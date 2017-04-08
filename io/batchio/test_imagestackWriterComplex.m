% Test the classes imagestackReaderComplex and imagestackWriterComplex.
%
% Generate a small number of complex-values matrices, save them using
% imagestackWriterComplex, read them back using imagestackReaderComplex and
% compare the original (in-memory) stack with the one that has been
% written-read. 
% The difference should be zero (both arrays should agree to the bit).
%
% Yoel Shkolnisky, April 2017.

% Generate projections.
n=65;
K=100;
images1=rand(n,n,K,'single')+1i.*ones(n,n,K,'single');

fname='tmp.mrc';
% Write to disk.
outstack=imagestackWriterComplex(fname,K,1,50);
outstack.append(images1);
outstack.close;

% Read saved images from MRC
imreader=imagestackReaderComplex(fname);

images2=zeros(n,n,K);
for idx=1:K
    images2(:,:,idx)=imreader.getImage(idx);
end

% Compare. projs and projs2 should agree to the bit.
err=norm(images1(:)-images2(:));
fprintf('err=%e\n',err);
if err==0
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end