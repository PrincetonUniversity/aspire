function [ masked ] = mask_corners( im )

% Tejal Bhamre, Oct 2015
%EMDB data not masked at corners, masking to zero.

L=size(im,1);
nim=size(im,3);
if mod(L,2)==1
    N=floor(L/2);
else 
    N=floor(L/2)-1;
end
[x, y]=meshgrid(-N:N, -N:N); 
r=sqrt(x.^2+y.^2);
r_max=floor(L/2)-1;

for i=1:nim
temp=zeros(L,L);
I=im(:,:,i);
temp(r<=r_max)=I(r<=r_max);
masked(:,:,i)=temp;
end

end

