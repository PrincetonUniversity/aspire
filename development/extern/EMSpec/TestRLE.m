% TestRLE
figure(1);
SetGrayscale;
nf=1;  % image complexity

n=256;
img=GaussFilt(randn(n),nf*1.33/n)>0;

% img=zeros(n,n);
% img(:,n/2+1:n)=1;
% img(:,n-3:n)=0;
% % img(:,n-2:n-1)=1;
% img(:,n-1:n)=1;
% 
% % img=1-img;

SetGrayscale;
subplot(121);
imacs(img);
code=BinaryRLEncode(img);
img2=BinaryRLDecode(code);
subplot(122);
imacs(img2);
title(numel(code.data));

% code.data
