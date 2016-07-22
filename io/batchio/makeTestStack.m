function makeTestStack(fname,sz,N)
%
% MAKETESTSTACK     Creare a stack of test images for test imagestack
%
% makeTestStack(fname,N)
%   Create a test image stack containing N images of size sz stored in the
%   MRC file fname.
%
% Example:
%   makeTestStack('test.mrc',65,1000);
%
% Yoel Shkolnisky, May 2014.

if isscalar(sz)
    sz=[sz sz];
end

stack=zeros(sz(1),sz(2),N);
for k=1:N
    im=makeTestImage(sz,k);
    stack(:,:,k)=im;
end
WriteMRC(stack,1,fname);
    
