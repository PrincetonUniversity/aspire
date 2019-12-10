
% Input:  2 sphere grids of equal size. 
% Output: A set of relative rotations Rij where Ri is from grid 1 and Rj
% from grid 2

function [Rijs]=genRijsFrom2GridsD2(grid1,grid2)
 %% Generate relative rotation quadruplets of D2 from a rotations grid
nrot=size(grid1,3);
npairs=nchoosek(nrot,2);
Rijs=zeros(3,3,4,npairs);

g_x=diag([1,-1,-1]);
g_y=diag([-1,1,-1]);
g_z=diag([-1,-1,1]);

idx=0;
for i=1:nrot-1
    for j=i+1:nrot
        idx=idx+1;
        Rijs(:,:,1,idx)=grid1(:,:,i)'*grid2(:,:,j);
        Rijs(:,:,2,idx)=grid1(:,:,i)'*g_x*grid2(:,:,j);
        Rijs(:,:,3,idx)=grid1(:,:,i)'*g_y*grid2(:,:,j);
        Rijs(:,:,4,idx)=grid1(:,:,i)'*g_z*grid2(:,:,j);
    end
end
   