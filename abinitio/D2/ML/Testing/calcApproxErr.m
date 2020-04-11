
function [sum_err,err_out,cls_err]=calcApproxErr(Rijs_gt,Rijs_est,cls_gt,cls_est)
%% Calculate error between approximated Rijs (approx) and Ground Truth (gt)
cls_err=[];
ps=perms(4:-1:1);
Rijs_est_perms=zeros([size(Rijs_est),24*2]);
Rijs_est_Jified=multi_Jify(Rijs_est);
for i=1:24
    Rijs_est_perms(:,:,:,:,i)=Rijs_est(:,:,ps(i,:),:);
    Rijs_est_perms(:,:,:,:,24+i)=Rijs_est_Jified(:,:,ps(i,:),:);
end
Rijs_dup=repmat(Rijs_gt,1,1,1,1,48);
s=size(Rijs_dup);
err=squeeze(reshape((Rijs_est_perms-Rijs_dup),[9,s(3:5)]));
err=squeeze(sqrt(sum((err.^2),1)));
[sum_err,best_var]=min(squeeze(sum(err,1)),[],2);
err_out=zeros(s(4),4);
for i=1:s(4)
    err_out(i,:)=err(:,i,best_var(i));
end

if nargin>2
    s=size(Rijs_est);
    cls_err=zeros(s(4),3);
    cls_est_perms=uint16(zeros([4,2,s(4),24*2]));
    cls_est_Jified=mod(cls_est+180,360);
    for i=1:24
        cls_est_perms(:,:,:,i)=cls_est(ps(i,:),:,:);
        cls_est_perms(:,:,:,24+i)=cls_est_Jified(ps(i,:),:,:);
    end
    cls_est_perms=double(cls_est_perms);
    
    cls_dup=double(repmat(cls_gt,1,1,1,48));
    s=size(cls_dup);
    ang_diff=squeeze(reshape(abs(cls_est_perms-cls_dup),[8,s(3:4)]));
    tmp=reshape(ang_diff,8*s(3)*s(4),1);
    more_than_180_idx=tmp>180;
    tmp(more_than_180_idx)=360-tmp(more_than_180_idx);
    ang_diff=reshape(tmp,[8,s(3:4)]);
    err_norm=squeeze(sqrt(sum((ang_diff.^2),1)));
    [err_norm,best_var_norm]=min(err_norm,[],2);
    err_max=squeeze(max(ang_diff,[],1));
    [err_max,best_var_max]=min(err_max,[],2);
    err_mean_abs=squeeze(sum(abs(ang_diff),1))/8;
    [err_mean_abs,best_var_mean_abs]=min(err_mean_abs,[],2);
    
    cls_err(:,1)=err_norm;
    cls_err(:,2)=err_max;
    cls_err(:,3)=err_mean_abs;
    
end

end