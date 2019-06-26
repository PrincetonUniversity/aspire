function [rot_alligned,err_in_degrees,mse] = analyze_results_Dn(rots,params)

if ~params.real_data && params.analyzeResults
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % step 13  : Results Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('\nStep 13: Analyzing results');
    [mse, rot_alligned, sign_g_Ri] = check_rotations_error(rots,params);
    fprintf('MSE of rotations estimate: %e\n',mse);
    
    err_in_degrees = check_degrees_error(rot_alligned,params);
    fprintf('median error in degrees: %e\n',median(err_in_degrees));
    fprintf('mean error in degrees: %e\n',mean(err_in_degrees));
    fprintf('std error in degrees: %e\n',std(err_in_degrees));
    fprintf('min error in degrees: %e\n',min(err_in_degrees));
    fprintf('max error in degrees: %e\n',max(err_in_degrees));
    
else
    rot_alligned = rots;
    mse = Inf;
    err_in_degrees = Inf;
end

end



function [mse, rot_alligned, sign_g_Ri] = check_rotations_error(rotations,params)
%
%
% explanation goes here...
%
% flag = 1 if no J-conjugacy, and 2 if there is a global J-conjugacy
%
% check the reconstruction error of R_1, ... ,R_K up to a global
% rotation/reflection O, and up to J-conjugate


% If reference rotations are given, compare the resluting rotations against
% true ones.

% 

% The estimation of every rotation is up to a rotation by pi. We don't know if we estimated R_i or gRi. 
%This ambiguity exist independently for every estimated rotation.
%Therefore, for error analysis we perform a global g-synchronization.
%Leaving only a possible,single, global rotation for all images. 
J=diag([1 1 -1]);
%global J;
g_s = zeros(3,3,params.FOUR);
g_s(:,:,1) = eye(3);
g_s(:,:,2) = diag([1 -1 -1]);
g_s(:,:,3) = diag([-1 1 -1]);
g_s(:,:,4) = diag([-1 -1 1]);

[sign_g_Ri] = g_sync(rotations,params);

    N=size(rotations,3);
    mse=zeros(2,6,8,2);
    %g=zeros(3,3,N,6,8);
    sign_mat=zeros(3,3,8);
    rot_alligned=zeros(3,3,N);
    diff=zeros(N,1);
    Rref=zeros(3,3,N);
    
    for k=1:N
        Rref(:,:,k)=q_to_rot(params.refq(:,k));
    end
    
    %generate all 3! row swap matrices
    row_swap_mat=zeros(3,3,6);
    idx=0;
    for i=1:3
        for j=1:i-1
            idx=idx+1;
            row_swap_mat(1,i,idx)=1;
            row_swap_mat(2,j,idx)=1;
            row_swap_mat(3,6-i-j,idx)=1;
        end
        for j=i+1:3
            idx=idx+1;
            row_swap_mat(1,i,idx)=1;
            row_swap_mat(2,j,idx)=1;
            row_swap_mat(3,6-i-j,idx)=1;
        end
    end
    
    %generate all 2^3 sign matrices
    idx=0;
    for i=[-1 1]
        for j=[-1 1]
            for k=[-1 1]
                idx=idx+1;
                sign_mat(:,:,idx)=diag([i j k]);
            end
        end
    end
    for n=1:2
        gs=sign_g_Ri(n,:);
        for alpha=[0 1]            
            for l=1:6
                for m=1:8        
                    for i=1:N
    %                     [g(:,:,i,l,m),d]=calcMinDist(row_swap_mat(:,:,l),...
    %                         sign_mat(:,:,m),Rref(:,:,i),rotations(:,:,i),alpha);
                        [d]=calcMinDist(row_swap_mat(:,:,l),...
                            sign_mat(:,:,m),gs(i),Rref(:,:,i),rotations(:,:,i),alpha);              
                        mse(n,l,m,alpha+1)=mse(n,l,m,alpha+1)+d;
                    end
                end
            end
        end
    end
    mse=mse/N;  
    
%     [M1,I]=min(mse(:,:,1),[],1);
%     [M1,m1]=min(M1);
%     l1=I(m1);    
%     [M2,I]=min(mse(:,:,2),[],1);
%     [M2,m2]=min(M2);
%     l2=I(m2);
%     if mse(l1,m1,1)<mse(l2,m2,2)
%         mse=M1;
%         alpha=0;
%         m=m1;
%         l=l1;
%     else
%         mse=M2;
%         alpha=1;
%         m=m2;
%         l=l2;
%     end
    
    [M,I]=min(mse(:));
    [n,l,m,alpha] = ind2sub(size(mse),I(1));
    alpha=alpha-1;
    mse=M;
    sign_g_Ri=sign_g_Ri(n,:);
    
    Ot=row_swap_mat(:,:,l);
    S=sign_mat(:,:,m);
    for i=1:N
        rot_alligned(:,:,i)=g_s(:,:,sign_g_Ri(i))*J^alpha*S*Ot*rotations(:,:,i)*J^alpha;
        diff(i)=norm(Rref(:,:,i)-rot_alligned(:,:,i),'fro');
    end        


% O = V*U.';
% rot_alligned = zeros(3,3,K);
% 
% % Plot estimation errors
% diff=zeros(K,1);
% for k=1:K
%     Rref=q_to_rot(params.refq(:,k)).';
%     
%     R= O*rotations(:,:,k);
%     if flag ==2
%         R = J*R*J;
%     end
%     rot_alligned(:,:,k) = R;
%     
%     diff(k) = norm(R - Rref,'fro');
% end
% mse= sum(diff.^2)/K;


return;

if params.debug
    figure,
    hist(diff, 40);
    set(gca,'FontSize', 20);
    title('hist of err_{Fro} of reconstructed R_i')    
end

end

% function [diff, mse,rot_alligned]=check_rotations_error(rotations,refq)
% 
%     N=size(rotations,3);
%     mse=zeros(6,8,2);
%     g=zeros(3,3,N,6,8);
%     sign_mat=zeros(3,3,8);
%     rot_alligned=zeros(3,3,N);
%     diff=zeros(N,1);
%     
%     %generate all 3! row swap matrices
%     row_swap_mat=zeros(3,3,6);
%     idx=0;
%     for i=1:3
%         for j=1:i-1
%             idx=idx+1;
%             row_swap_mat(1,i,idx)=1;
%             row_swap_mat(2,j,idx)=1;
%             row_swap_mat(3,6-i-j,idx)=1;
%         end
%         for j=i+1:3
%             idx=idx+1;
%             row_swap_mat(1,i,idx)=1;
%             row_swap_mat(2,j,idx)=1;
%             row_swap_mat(3,6-i-j,idx)=1;
%         end
%     end
%     
%     %generate all 2^3 sign matrices
%     idx=0;
%     for i=[-1 1]
%         for j=[-1 1]
%             for k=[-1 1]
%                 idx=idx+1;
%                 sign_mat(:,:,idx)=diag([i j k]);
%             end
%         end
%     end
%                     
%     for alpha=[0 1]            
%         for l=1:6
%             for m=1:8        
%                 for i=1:N
%                     [g(:,:,i,l,m),d]=calcMinDist(row_swap_mat(:,:,l),...
%                         sign_mat(:,:,m),refq(:,:,i),rotations(:,:,i),alpha);
%                     mse(l,m,alpha+1)=mse(l,m,alpha+1)+d;
%                 end
%             end
%         end
%     end
%     mse=mse/N;  
%     
%     [M1,I]=min(mse(:,:,1),[],1);
%     [M1,m1]=min(M1);
%     l1=I(m1);    
%     [M2,I]=min(mse(:,:,2),[],1);
%     [M2,m2]=min(M2);
%     l2=I(m2);
%     if mse(l1,m1,1)<mse(l2,m2,2)
%         mse=M1;
%         alpha=0;
%         m=m1;
%         l=l1;
%     else
%         mse=M2;
%         alpha=1;
%         m=m2;
%         l=l2;
%     end
%     
%     Ot=row_swap_mat(:,:,l);
%     S=sign_mat(:,:,m);
%     J=diag([1 1 -1]);
%     for i=1:N
%         rot_alligned(:,:,i)=J^alpha*S*Ot*rotations(:,:,i)*J^alpha;
%         diff(i)=norm(refq(:,:,i)-rot_alligned(:,:,i),'fro');
%     end        
% end

%check the best 
function [mindist]=calcMinDist(row_swap_mat,S,sign_g_Ri,Ri,Ri_est,alpha)
%function [g,mindist]=calcMinDist(row_swap_mat,S,Ri,Ri_est,alpha)    
    J=diag([1 1 -1]);
    %mindist=inf;
    g_alpha=[1 1 1;1 -1 -1;-1 1 -1;-1 -1 1];
    
    d=norm(Ri-diag(g_alpha(sign_g_Ri,:))*J^alpha*S*...
        row_swap_mat*Ri_est*J^alpha,'fro');

    mindist=d^2;
end

function err_in_degrees = check_degrees_error(rots,params)

%global g;
nImages = size(params.refq,2);
assert(size(rots,3) == nImages);
n_theta = params.n_theta;
d_theta = 2*pi/params.n_theta;

% precalculate g^s, s=0..params.FOUR-1
g_s = zeros(3,3,params.FOUR);
g_s(:,:,1) = eye(3);
g_s(:,:,2) = diag([1 -1 -1]);
g_s(:,:,3) = diag([-1 1 -1]);
g_s(:,:,4) = diag([-1 -1 1]);

for i=1:nImages
    rots(:,:,i)=rots(:,:,i)';
end

err_in_degrees = zeros(1,n_theta*nImages);
for k = 1:nImages
    rays_gt = zeros(n_theta,3);
    Rk_gt  = q_to_rot(params.refq(:,k))';
    for j = 1:n_theta
        rays_gt(j,:) = cos((j-1)*d_theta).*Rk_gt(:,1) ...
            +  ...
            sin((j-1)*d_theta).*Rk_gt(:,2);
    end
    
    c_err_in_degrees = cell(1,4);
    for s = 1:params.FOUR
        rays = zeros(n_theta,3);
        Rk   = g_s(:,:,s)*rots(:,:,k);
        for j = 1:n_theta
            rays(j,:) = cos((j-1)*d_theta).*Rk(:,1) ...
                +  ...
                sin((j-1)*d_theta).*Rk(:,2);
        end
        
        cos_angle   = rays(:,1).*rays_gt(:,1) + rays(:,2).*rays_gt(:,2) + rays(:,3).*rays_gt(:,3);
        cos_angle(cos_angle>1)=1;
        cos_angle(cos_angle<-1)=-1;
        err         = acos(cos_angle);
        c_err_in_degrees{s} = err * 180/pi;
    end
    meanErrs   = cellfun(@mean, c_err_in_degrees);
    [~, opt_s] = min(meanErrs);
    err_in_degrees((k-1)*n_theta+1:k*n_theta) = c_err_in_degrees{opt_s};
end


% err=acos(cos_angle);
hist(err_in_degrees,50);
h=gcf;
set(gca,'FontSize',14);
%title('Estimation error (in degrees)','FontSize',16);


end

function [classes] = g_sync(rotations, params)
%function [rotations_g_synced,sign_g_Ri] = g_sync(rotations, params)

% The in-plain rotation stage retrieves all rotation matrices. However, every calculated rotation might 
% be a the rotated version of g^s, s=1,...,n from the ground truth. 
%Moreover, the rotation to be applied might be differernt for every image.
% We therefore synchronize all rotation so that only a possibly single global
% rotation should be applied to all rotation.
% . Since the molecule is symmetric this doesn't matter as far as the reconstruction is concerned, but it is important if 
% wish to compare to ground truth

nImages = size(rotations,3);

J=diag([1 1 -1]);
%global g;

gs = zeros(3,3,4);
gs(:,:,1) = eye(3);
gs(:,:,2) = diag([1 -1 -1]);
gs(:,:,3) = diag([-1 1 -1]);
gs(:,:,4) = diag([-1 -1 1]);

Rs_t = zeros(3,3,nImages);
for i = 1:nImages
    Rs_t(:,:,i) = q_to_rot(params.refq(:,i));
end

% J_Rs_J_t = zeros(3,3,nImages);
% for i = 1:nImages
%     J_Rs_J_t(:,:,i) = J*Rs_t(:,:,i)*J;
% end

J_Rs_J= zeros(3,3,nImages);
for i = 1:nImages
    J_Rs_J(:,:,i) = J*rotations(:,:,i)*J;
end

A_g = zeros(2*nImages,2*nImages);
%A_g2 = zeros(nImages,nImages);

for i = 1:nImages-1
    

    %RiRjs_t = multiprod(Ri_t.', Rjs_t);
    Ri    = rotations(:,:,i);
    Rjs   = rotations(:,:,i+1:end);
    RiRjs = multiprod(Ri.', Rjs);   
    J_Ri_J    = J_Rs_J(:,:,i);
    J_Rjs_J   = J_Rs_J(:,:,i+1:end);
    J_RiRjs_J = multiprod(J_Ri_J.', J_Rjs_J);   
    
%   RiRjs  = zeros(3,3,nImages-i,params.FOUR);
%     for s  = 1:params.FOUR
%         g_s = gs(:,:,s);
%         RiRjs(:,:,:,s) = multiprod(Ri.'*g_s, Rjs);
%     end
    Ri_t      = Rs_t(:,:,i);
    Rjs_t    = Rs_t(:,:,i+1:end);
    RiRjs_t  = zeros(3,3,nImages-i,params.FOUR);
    for s  = 1:params.FOUR
        g_s = gs(:,:,s);
        RiRjs_t(:,:,:,s) = multiprod(Ri_t.'*g_s, Rjs_t);
    end
    
    % compute the squared frobenious norm
    norm_diff      = (bsxfun(@minus,RiRjs,RiRjs_t)).^2;
    J_norm_diff_J = (bsxfun(@minus,J_RiRjs_J,RiRjs_t)).^2;
    
    norm_diff     = reshape(norm_diff,9,nImages-i,params.FOUR);
    J_norm_diff_J = reshape(J_norm_diff_J,9,nImages-i,params.FOUR);
    
    % we don't care if there is conjugation or not - all we want is to know
    % s for g^s
    norm_diff_concat = cat(3,norm_diff,J_norm_diff_J);
    
    [diff,ii] = min(sum(norm_diff_concat),[],3);
    
    % remove the effect of concatination above.
    ii = mod(ii,4);
    
    rep=zeros(2,2*length(ii));
    idx=find(ii==0);
    rep(1,2*idx-1)=-1;
    rep(2,2*idx)=-1;
    
    idx=find(ii==1);
    rep(1,2*idx-1)=1;
    rep(2,2*idx)=1;

    idx=find(ii==2);
    rep(1,2*idx-1)=1;
    rep(2,2*idx)=-1;

    idx=find(ii==3);
    rep(1,2*idx-1)=-1;
    rep(2,2*idx)=1;


    
    A_g(2*i-1:2*i,2*i+1:end) = rep;
    
    
    
%     tmp1=ii;
%     tmp2=ii;
%     tmp1(tmp1<2)=1; % 0,1 ->1
%     tmp1(tmp1>1)=-1; % 2,3 -> -1
%     tmp2(tmp2==1)=-1; % 1 -> -1
%     tmp2(tmp2==2)=-1; % 2 -> -1
%     tmp2(mod(tmp2,3)==0)=1; %3-> 1, 0->1
%     A_g(i,i+1:end) = tmp1;%exp(-sqrt(-1)*params.alpha*(ii-1)); 
%     A_g2(i,i+1:end)=tmp2;
    
end
% A_g(k,l) is exp(-j(-theta_k+theta_j)) so we use transpose and not
% conjugate-transpose to obtain lower triangular entries
A_g = A_g + A_g';
%A_g2 = A_g2 + A_g2';
% Diagonal elements correspond to exp(-i*0) so put 1. 
%This is important only for debug/verification purposes that spectrum is (K,0,0,0...,0)
A_g = A_g + eye(size(A_g)); 
%A_g2 = A_g2 + eye(size(A_g2)); 

[v,d] = eigs(A_g,10);% eigs(A_g, 10, 'lm');
[evals, ind] = sort(diag(abs(d)), 'descend');
evecs = v(:,ind(1:2));
% [v2,d2] = eig(A_g2);% eigs(A_g, 10, 'lm');
% [evals2, ind2] = sort(diag(d2), 'ascend');
% evect2 = v2(:,ind2(1));

evect1=sign(evecs(1:2:end-1,1));
evect2=sign(evecs(2:2:end,2));
classes=zeros(2,nImages);
classes(1,evect1==1 & evect2==1)=1;
classes(1,evect1==1 & evect2==-1)=2;
classes(1,evect1==-1 & evect2==-1)=4;
classes(1,evect1==-1 & evect2==1)=3;

classes(2,evect1==1 & evect2==1)=1;
classes(2,evect1==1 & evect2==-1)=2;
classes(2,evect1==-1 & evect2==-1)=4;
classes(2,evect1==-1 & evect2==1)=3;

% c1=evect1;
% c1(evect1<0)=2;
% idxs=find(evect1>0);
% tmp=evect1(idxs);
% tmp(evect2(idxs)<0)=4;
% c1(idxs)=tmp;
% idxs=find(evect1<0);
% tmp=evect1(idxs);
% tmp(evect2(idxs)<0)=3;
% tmp(evect2(idxs)>0)=2;
% c1(idxs)=tmp;
% classes=c1';
% 
% classes=zeros(8,length(evect1));
% classes(1,:)=c1;
% classes(2,:)=c1;
% classes(3,:)=c1;
% classes(4,:)=c1;
% classes(2,find((c1==1)+(c1==4)))=-c1(find((c1==1)+(c1==4)))+5;
% classes(3,find((c1==2)+(c1==3)))=-c1(find((c1==2)+(c1==3)))+5;
% classes(4,:)=-c1+5;
% 
% idxs1=find((c1==1)+(c1==4));
% idxs2=find((c1==2)+(c1==3));
% c2=c1;
% c2(idxs1)=(c1(idxs1)+5)/3;
% c2(idxs2)=3*c1(idxs2)-5;
% classes(5,:)=c2;
% classes(6,:)=c2;
% classes(7,:)=c2;
% classes(8,:)=c2;
% classes(6,find((c1==1)+(c1==4)))=-c1(find((c1==1)+(c1==4)))+5;
% classes(7,find((c1==2)+(c1==3)))=-c1(find((c1==2)+(c1==3)))+5;
% classes(8,:)=-c1+5;

disp(['Outer g syn. First 5 eigenvalues : ' num2str(evals(1:4)')]);

% angles = exp(sqrt(-1)*params.alpha*(0:params.FOUR-1))';
% rotations_g_synced = zeros(3,3,nImages);
%sign_g_Ri = zeros(nImages,1);
% for ii  = 1:nImages
%     zi  = evect1(ii);
%     zi  = zi/abs(zi); % rescale so it lies on unit circle
%     signs(ii)=zi;
%     
%     % matlab returnes complex number phases (angles) in [-pi,+pi]. Since a ccw and a cw closest are just as good, 
%     % we take the absolute value of the angle
%     angleDists = abs(angle(zi./angles));
%     
%     [~,I] = min(angleDists);
%     sign_g_Ri(ii) = params.FOUR-(I(1)-1);
% %     R = rotations(:,:,ii);
% %     rotations_g_synced(:,:,ii) = gs(:,:,sign_g_Ri(ii))*R;    
% end

end
