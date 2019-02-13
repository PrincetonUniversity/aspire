
%Compute the group members of Dn
function [g1,g2]=gns(n)
g1=zeros(3,3,n);
g2=zeros(3,3,n);
theta=2*pi/n;
g_z=[cos(theta),-sin(theta),0;
     sin(theta),cos(theta),0;
     0,0,1];
for j=0:n-1
   g1(:,:,j+1)=g_z^j; 
   g2(:,:,j+1)=axisang2rot([cos(pi*j/n);sin(pi*j/n);0],pi);
end
 
