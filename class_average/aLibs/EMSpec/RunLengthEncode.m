function code=RunLengthEncode(img)
% function code=RunLengthEncode(img)
% Use run-length coding to encode the binary image img.  Run lengths are
% encoded by 1 or more bytes (2^7 per byte).  This scheme is good for
% sparse images. The structure code is decoded by the function BinaryRLD.
% F. Sigworth 30 Jun 12

code.sz=size(img);
code.version=1;
z=img(:);
dz=[z(1); diff(z)];
runs=find(dz~=0);
nr=numel(runs);
if nr==0
    code.data=[0 0];
    return
end;

druns=[runs(1)-1; diff(runs(1:nr))];
runenc=uint8(4*nr);
k=1;
for j=1:nr
    [bytes nb]=EncodeRun(druns(j));
    runenc(k:k+nb-1)=bytes(1:nb);
    k=k+nb;
end;
code.data=runenc(1:k-1);


function [bytes nb]=EncodeRun(r)
bytes=uint8(8);
r1=r;
nb=1;
for i=1:8
    if r1<128
        bytes(nb)=r1;
        r1=0;
        break;
    else
        bytes(nb)=128+mod(r1,128);
        r1=floor(r1/128);
        nb=nb+1;
    end;
end;
if r1>0
    error(['Overflow in encoding ' num2str(r)]);
end;
