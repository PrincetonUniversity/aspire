% function verifyImage(im);
%
% Verify input image.
%
% The function verifies that the input image is a 3D image of size nxnxn
% with n even. Otherwise the function aborts with an error message.
%
% Yoel Shkolnisky

function verifyImage(im);

% Verify that the input is a 3D image of size nxnxn
s=size(im);
if length(s) ~= 3
   error('Input must be a 3D image');
end

if (s(1)-s(2)~=0) | (s(2)-s(3)~=0)
   error('Input image must be cube');
end

if (mod(s(1),2)~=0) | (mod(s(2),2)~=0) | (mod(s(3),2)~=0)
   error('Input image must have even sides');
end
