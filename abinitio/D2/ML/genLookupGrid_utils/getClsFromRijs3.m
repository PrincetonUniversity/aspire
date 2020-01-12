
%% Input:  cls = pairs of indices (in degrees) of 4-tuples of common lines
%  Output: common lines
function [cls_idx,sub_idx]=getClsFromRijs3(cls)
%% Convert to linear indexs corresponding to correlation tables between
%  image rays in Fourier space. 
rowSub=mod(round(cls(:,1)),360);
colSub=mod(round(cls(:,2)),360);

% Restrict Rj in-plane coordinates to <180 degrees
is_geq_than_pi=(colSub>=180);
if sum(is_geq_than_pi)>0
    colSub(is_geq_than_pi)=colSub(is_geq_than_pi)-180;
    rowSub(is_geq_than_pi)=bsxfun(@mod,rowSub(is_geq_than_pi)+180,360);
end

%Add +1 since indices in MATLAB start from 1 (i degrees ray has i-1 index).
rowSub=rowSub+1;
colSub=colSub+1;
sub_idx=[rowSub,colSub];

%convert to linear indexes in 360*180 correlation matrix
cls_idx=uint16(sub2ind([360,180],rowSub,colSub));

end