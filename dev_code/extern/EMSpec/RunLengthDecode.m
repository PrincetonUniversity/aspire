function img=RunLengthDecode(code)
% function img=RunLengthDecode(code)
% Perform run-length decoding to recreate a binary image.
% Decodes the data structure created by BinaryRLE().
% F. Sigworth 30 Jun 12
% if code.version ~=1
%     error(['Unrecognized RLE version: ' num2str(code.version)]);
% end;
img=single(zeros(code.sz));
z=img(:);
n=numel(z);
nc=numel(code.data);
k=1;
j=1;
v=0;
while k<=nc
    [r k]=DecodeRun(code.data,k);
%     disp([k r j])
%     r
    img(j:j+r-1)=v;
    v=1-v;
    j=j+double(r);
%     j
end;
img(j:n)=v;
end


function [r k]=DecodeRun(z,k)
r=0;
m=1;
r1=double(z(k));
while r1>127
    r1=r1-128;
    r=r+m*r1;
    k=k+1;
    m=m*128;
    r1=double(z(k));
end;
r=r+r1*m;
k=k+1;
end
