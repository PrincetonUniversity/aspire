
function [plane_proj]=projectToAllPlanes(grid)

xy_proj=grid;
xy_proj(3,:)=0;
norms=sqrt(sum(xy_proj.^2,1));
xy_proj=xy_proj./repmat(norms,3,1);

xz_proj=grid;
xz_proj(2,:)=0;
norms=sqrt(sum(xz_proj.^2,1));
xz_proj=xz_proj./repmat(norms,3,1);

yz_proj=grid;
yz_proj(1,:)=0;
norms=sqrt(sum(yz_proj.^2,1));
yz_proj=yz_proj./repmat(norms,3,1);
plane_proj=cat(3,xy_proj,xz_proj,yz_proj);

