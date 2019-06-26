
%Input:  Relative rotations Rijs
%Output: common lines
function [cls_idx,sub_idx]=getClsFromRijs_scls_Dn(cls)

%Convert to linear indexs corresponding to correlations lookup table
%generated in cryo_clmatrix
rowSub=mod(cls(:,1),360);
colSub=mod(cls(:,2),360);

%Add +1 since indexes start from 1 (i degrees ray has i-1 index).
rowSub=rowSub+1;
colSub=colSub+1;

% Restrict Rj in-plane coordinates to <180 degrees
is_geq_than_pi=(colSub>180);
if sum(is_geq_than_pi)>0
    colSub(is_geq_than_pi)=colSub(is_geq_than_pi)-180;
    rowSub(is_geq_than_pi)=bsxfun(@mod,rowSub(is_geq_than_pi)+180,360);
end
rowSub(rowSub==0)=360;

%Add +1 since indexes start from 1 (i degrees ray has i-1 index).
%rowSub=rowSub+1;
%colSub=colSub+1;
sub_idx=[rowSub,colSub];

%converst to linear indexes in 360*180 correlation matrix
cls_idx=uint16(sub2ind([360,180],rowSub,colSub));

end