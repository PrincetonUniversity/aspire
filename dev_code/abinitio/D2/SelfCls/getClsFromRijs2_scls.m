
%Input:  Relative rotations Rijs
%Output: common lines
function [cls_idx,sub_idx]=getClsFromRijs2_scls(cls)

%npairs=size(cls_in,3);
% cls(:,1)=atan2(Rijs(1,3,:),-Rijs(2,3,:));
% cls(:,2)=atan2(-Rijs(3,1,:),Rijs(3,2,:));

% cls=mod(cls+2*pi,2*pi); %make all angles non negative 
% cls=cls*180/pi; %convert to degrees 

%Convert to linear indexs corresponding to correlations lookup table
%generated in cryo_clmatrix
rowSub=mod(cls(:,1)-1,360);
colSub=mod(cls(:,2)-1,360);

% Restrict Rj in-plane coordinates to <180 degrees
is_geq_than_pi=(colSub>=180);
if sum(is_geq_than_pi)>0
    colSub(is_geq_than_pi)=colSub(is_geq_than_pi)-180;
    rowSub(is_geq_than_pi)=bsxfun(@mod,rowSub(is_geq_than_pi)+180,360);
end

%At this point all thetas of Rj are >=0 and <180. We cannot just round
%up j thetas >179.5 since 0 angle is not equivalent to 180. so we must
%floor it. For the respective i's we subtract the change in j angle and
%then round to minimize that error in i angle to be <=0.5 degree

% Old Code
% colSub_close_to_180_idx=colSub>179.5;
% diff_after_floor=colSub(colSub_close_to_180_idx)-...
%     floor(colSub(colSub_close_to_180_idx));
% colSub(colSub_close_to_180_idx)=floor(colSub(colSub_close_to_180_idx));
% rowSub(colSub_close_to_180_idx)=...
%     round(rowSub(colSub_close_to_180_idx)-diff_after_floor);
% rowSub=round(rowSub);
% colSub=round(colSub);
% rowSub=mod(rowSub+360,360);

% New code
% colSub_close_to_180_idx=colSub>179.5;
% colSub(colSub_close_to_180_idx)=0;
% rowSub(colSub_close_to_180_idx)=mod(rowSub(colSub_close_to_180_idx)+180,360);
% rowSub=round(rowSub);
% colSub=round(colSub);
% rowSub=mod(rowSub+360,360);

%Add +1 since indexes start from 1 (i degrees ray has i-1 index).
rowSub=rowSub+1;
colSub=colSub+1;
sub_idx=[rowSub,colSub];

%converst to linear indexes in 360*180 correlation matrix
cls_idx=uint16(sub2ind([360,180],rowSub,colSub));

end