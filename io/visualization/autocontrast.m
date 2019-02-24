function im2=autocontrast(im)

s=size(im);
im=im(:);
frac=0.95; %remove the 5% top values and 5% lowest values;
rf=(1-frac)/2; % fraction of values to chop from each side of the histogram

[~,I]=sort(im);
imlen=prod(s);

Ilow=ceil(imlen*rf);
Ihigh=imlen-Ilow;
im(I(1:Ilow))=im(I(Ilow));
im(I(Ihigh:imlen))=im(I(Ihigh));
im2=reshape(im,s(1),s(2));
