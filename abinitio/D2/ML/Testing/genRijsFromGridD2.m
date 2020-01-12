
function [Rijs,cls]=genRijsFromGridD2(grid)
%% Generate relative rotation quadruplets of D2 from a rotations grid
nrot=size(grid,3);
npairs=nchoosek(nrot,2);
Rijs=zeros(3,3,4,npairs);
cls=zeros(4,2,npairs);

g_x=diag([1,-1,-1]);
g_y=diag([-1,1,-1]);
g_z=diag([-1,-1,1]);

idx=0;
for i=1:nrot-1
    for j=i+1:nrot
        idx=idx+1;
        Rijs(:,:,1,idx)=grid(:,:,i)'*grid(:,:,j);
        Rijs(:,:,2,idx)=grid(:,:,i)'*g_x*grid(:,:,j);
        Rijs(:,:,3,idx)=grid(:,:,i)'*g_y*grid(:,:,j);
        Rijs(:,:,4,idx)=grid(:,:,i)'*g_z*grid(:,:,j);
        for l=1:4
            cls(l,1,idx)=mod(atan2(Rijs(1,3,l,idx),-Rijs(2,3,l,idx))+2*pi,2*pi);
            cls(l,2,idx)=mod(atan2(-Rijs(3,1,l,idx),Rijs(3,2,l,idx))+2*pi,2*pi);
        end
    end
end
cls=round(180*cls/pi);
cls1=squeeze(cls(:,1,:));
cls1=cls1(:);
cls1(cls1==0)=360;
cls2=squeeze(cls(:,2,:));
cls2=cls2(:);
cls2(cls2==0)=360;
geq_than_180_idx=cls2>180;
cls2(geq_than_180_idx)=cls2(geq_than_180_idx)-180;
cls1(geq_than_180_idx)=mod(cls1(geq_than_180_idx)+180,360);
cls1(cls1==0)=360;
cls(:,1,:)=reshape(cls1,4,npairs);
cls(:,2,:)=reshape(cls2,4,npairs);






