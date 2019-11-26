function im2=autocontrast(im)

if ndims(im)==3
    N=size(im,3);
else
    N=1;
end

im2=zeros(size(im));

for k=1:N
    tmp=im(:,:,k);
    s=size(tmp);
    
    frac=0.9; %remove the 5% top values and 5% lowest values;
    rf=(1-frac)/2; % fraction of values to chop from each side of the histogram
    
    [~,I]=sort(tmp);
    imlen=prod(s);
    
    Ilow=ceil(imlen*rf);
    Ihigh=imlen-Ilow;
    tmp(I(1:Ilow))=tmp(I(Ilow));
    tmp(I(Ihigh:imlen))=tmp(I(Ihigh));
    im2(:,:,k)=reshape(tmp,s(1),s(2));
end