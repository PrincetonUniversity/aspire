
% Filter viewing directions to leave only the slice of 180/n of the sphere
function oct=filterGridDn(grid,n)

xyProj=grid(1:2,:);
% norms=sqrt(sum(xyProj.^2,1));
% xyProj=xyProj./norms;
maxTheta=pi/n;
dirArgs=atan2(xyProj(2,:)',xyProj(1,:)');
oct=grid(:,(dirArgs<=maxTheta) & (dirArgs>=0));

if mod(n,2)>0
    theta=-pi/(2*n);
     Rz_theta=[cos(theta),-sin(theta),0;
     sin(theta),cos(theta),0;
     0,0,1];
    oct=Rz_theta*oct;
end




