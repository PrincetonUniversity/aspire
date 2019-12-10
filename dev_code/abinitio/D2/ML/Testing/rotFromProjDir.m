
%Input:  projDir = projection direction
%        inplane_ang = angle for in plane rotation about projDir
%Output: rotation matrix with projDir as third column

function R=rotFromProjDir(projDir,inplane_ang)
v=projDir;
v=v/norm(v);
a=inplane_ang;
u=[-v(2);v(1);0];
R=[cross(u,v)/norm(cross(u,v)),u,v];
Rz=[cos(a),-sin(a),0;
    sin(a),cos(a),0;
    0,0,1];
R=R*Rz;