
%Input:  If isempty(Q) then R is all 4 relative rotations R^Tg_kQ
%Output: cls = the coordinates of the common lines induced by the realtive 
%        rotations in degrees. cryo_plot_cls produces a plot of these
%        common lines. 
function [cls,corrs,shifts,Ax1,Ax2]=cryo_plot_cls2(R,Q,C,S)
cls=zeros(4,2);
cls_deb=zeros(4,2);
if nargin>1
    gs=cat(3,eye(3),diag([1,-1,-1]),diag([-1,1,-1]),diag([-1,-1,1]));
    for k=1:4        
        RQ=R'*gs(:,:,k)*Q;
        f=1;%sqrt(1-RQ(3,3)^2);
        cls(k,2)=atan2(-RQ(3,1,:,:)/f,RQ(3,2,:,:)/f);
        cls(k,1)=atan2(RQ(1,3,:,:)/f,-RQ(2,3,:,:)/f);
        
        %DEBUG
        q=cross(R(:,3),gs(:,:,k)*Q(:,3));
        q=q/norm(q);
        c2=(gs(:,:,k)*Q)'*q;
        c1=R'*q;
        cls_deb(k,2)=atan2(c2(2),c2(1));
        cls_deb(k,1)=atan2(c1(2),c1(1));
        
    end
else 
    % R is actually all R^Tg_kQ
    for k=1:4
        RQ=R(:,:,k);
        cls(k,2)=atan2(-RQ(3,1,:,:),RQ(3,2,:,:));
        cls(k,1)=atan2(RQ(1,3,:,:),-RQ(2,3,:,:));
    end
end
cls=mod(cls+2*pi,2*pi); %make all angles non negative 
cls_in_deg=round(cls*180/pi);
% bigger_than_180_idx=cls_in_deg(:,1)>=180;
% cls_in_deg(bigger_than_180_idx,1)=cls_in_deg(bigger_than_180_idx,1)-180;
% cls_in_deg(bigger_than_180_idx,2)=mod(cls_in_deg(bigger_than_180_idx,2)+180,360);
% cls(bigger_than_180_idx,1)=cls(bigger_than_180_idx,1)-pi;
% cls(bigger_than_180_idx,2)=mod(cls(bigger_than_180_idx,2)+pi,2*pi);
[eq_idx,~]=cryo_detect_equator_rots(R,5);
isscl=eq_idx>0 & norm(R-Q)<1e-10;
 start_idx=1;
if isscl
    start_idx=2;
end

figure;
subplot(1,2,1);
colors='brgy';
cls_labels=cell(4,1);
for k=start_idx:4    
    theta=[cls(k,1),cls(k,1)];
    rho=[0,1];    
    if (k>1 && ~isscl) || (k>2)
        hold on
    end
    polarplot(theta,rho,colors(k));
    cls_labels{k}=num2str(cls_in_deg(k,1));
end
hold off
Ax=gca;
Ax.RTickLabel=[];
Ax.RGrid='off';
Ax.ThetaGrid='off';
%[sorted_ang,ang_ord]=sort(round(cls(:,1)*180/pi));
[sorted_ang,ang_ord]=unique(round(cls(:,1)*180/pi));
thetaticks(sorted_ang);
thetaticklabels(cls_labels(ang_ord));
Ax1=Ax;

subplot(1,2,2);
for k=start_idx:4    
    theta=[cls(k,2),cls(k,2)];
    rho=[0,1];    
    if (k>1 && ~isscl) || (k>2)
        hold on
    end
    polarplot(theta,rho,colors(k));
    cls_labels{k}=num2str(cls_in_deg(k,2));
end
hold off
Ax=gca;
Ax.RTickLabel=[];
Ax.RGrid='off';
Ax.ThetaGrid='off';
%[sorted_ang,ang_ord]=sort(cls_in_deg(:,2));
[sorted_ang,ang_ord]=unique(cls_in_deg(:,2));
thetaticks(sorted_ang);
thetaticklabels(cls_labels(ang_ord));
Ax2=Ax;

cls=cls*180/pi;
cls=mod(round(cls),360)+1;
geq_180_idx=cls(:,2)>180;
cls(geq_180_idx,1)=mod(cls(geq_180_idx,1)+180,360);
cls(geq_180_idx,2)=cls(geq_180_idx,2)-180;

%If input contains correlations matrix output cls correlations
corrs=zeros(4,1);
if nargin>2    
    corrs(1)=C(cls(1,1),cls(1,2));
    corrs(2)=C(cls(2,1),cls(2,2));
    corrs(3)=C(cls(3,1),cls(3,2));
    corrs(4)=C(cls(4,1),cls(4,2));
end


%If input contains 1D shifts matrix output cls shifts
shifts=zeros(4,1);
if nargin>3
    shifts(1)=S(cls(1,1),cls(1,2));
    shifts(2)=S(cls(2,1),cls(2,2));
    shifts(3)=S(cls(3,1),cls(3,2));
    shifts(4)=S(cls(4,1),cls(4,2));
end
%cls=round(cls*180/pi); %convert to degrees 
%fig=gcf;





