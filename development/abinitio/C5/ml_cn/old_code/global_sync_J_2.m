function [vijs_out,viis_out,sign_ij_J,sign_ii_J] = global_sync_J_2(refq)
%
% Global J-synchronization of all third-row outer products. Specifically, the input 3n-by-3n matrix V whose
% (i,j)-th block (1<=i<j<=n) of size 3-by-3 is the estimate Rij for
% 1) Either $Ri^{T}g^{sij}Rj$, or $J Ri^{T}g^{sij}Rj J$ and sij is either 0,1,2 or 3 (if i!=j),
% 2) Either $Ri^{T}g^{si}Rj$, or $J Ri^{T}g^{si}Rj J$   and si is either      1 or 3 (if i==j)
% The method returns a matrix of same size with either all estimates have a
% spurious J in them, or none do at all.
%
% Input parameters:
%   vijs           A 3x3xn_choose_2 array where each slice holds an estimate for
%                  the corresponding outer-product vi*vj^{T} between the
%                  third rows of matrices Ri and Rj. Each such estimate
%                  might have a spurious J independently of other estimates
%   viis           A 3x3xn array where the i-th slice holds an estimate for
%                  the outer-product vi*vi^{T} between the
%                  third row of matrix Ri with itself. Each such estimate
%                  might have a spurious J independently of other estimates
% Output parameters:
%   vijs_out     A 3x3xn_choose_2 array where each slice holds an estimate for
%                the corresponding outer-product vi*vj^{T} 
%   sign_ij_J    An array of length n-choose-2 where each entry equals 1
%                or -1 depending whether the corresponding input 
%                estimate vivj had a spurious J or not.
%   sign_ii_J    An array of length n where each the i-thentry equals 1
%                or -1 depending whether the corresponding input 
%                estimate vii had a spurious J or not.


log_message('Global J synchronization');
% partition all off-diagonal blocks into two classes: Those with a spurious
% J in them and those without a spurious J in them.
nImages = size(refq,2);
% assert(size(vijs,3) == nchoosek(nImages,2));

if true
    viis = zeros(3,3,nImages);
    for i=1:nImages
        Ri_gt = q_to_rot(refq(:,i)).';
        viis(:,:,i) = Ri_gt(3,:).'*Ri_gt(3,:);
    end
end

J = diag([1,1,-1]);
vijs_new_cell = cell(nImages,nImages);
for i=1:nImages
    for j=i+1:nImages
        for k=j+1:nImages
            
            score_rank1 = zeros(1,32);
            s = 0;
            for si = 0:1
                for sj = 0:1
                    for sk = 0:1
                                    
                        vii = viis(:,:,i);
                        vjj = viis(:,:,j);
                        vkk = viis(:,:,k);
                        
                        vii = J^si*vii*J^si;
                        vjj = J^sj*vjj*J^sj;
                        vkk = J^sk*vkk*J^sk;
                        
                        
                        vij = vii*vjj;
                        vjk = vjj*vkk;
                        vki = vkk*vii;
                        
                        for sij=0:1
                            for sjk=0:1
                                
                                s = s + 1;
                                
                                vij = J^(sij)*vij*J^(sij); vij = vij/norm(vij,'fro');
                                vjk = J^(sjk)*vjk*J^(sjk); vjk = vjk/norm(vjk,'fro');
                                
                                V = [vii  vij vki';...
                                    vij' vjj vjk;...
                                    vki vjk' vkk];
                                
                                svals = svd(V/3);
                                
                                score_rank1(s) = norm(svals-[1,zeros(1,8)]',2);
                                
                            end
                        end
                    end
                end
            end
            [~,II] = min(score_rank1);
            s = num2str(dec2bin(II,5));
            
            if s(1) == 1
                vii = J*viis(:,:,i)*J;
            else
                vii = viis(:,:,i);
            end
            
            if s(2) == 1
                vjj = J*viis(:,:,j)*J;
            else
                vjj = viis(:,:,j);
            end
            
            if s(3) == 1
                vkk = J*viis(:,:,k)*J;
            else
                vkk = viis(:,:,k);
            end
            
            if s(4) == 1
                vij = J*vii*vjj*J;
            else
                vij = vii*vjj;
            end
            
            if s(5) == 1
                vjk = J*vjj*vkk*J;
            else
                vjk = vjj*vkk;
            end
            
            if s(6) == 1
                vki = J*vkk*vii*J;
            else
                vki = vkk*vii;
            end
        
            vijs_new_cell{i,j}(end+1) = vij;
            vijs_new_cell{j,k}(end+1) = vjk;
            vijs_new_cell{i,k}(end+1) = vki';
        end
    end
end


% next we need to choose a rep for each i,j. cluster to different groups
% and average over one of the groups (doesn't matter which one we choose)
vijs_averaged = zeros(3,3,nchoosek(nImages,2));
ind =0;
for i=1:nImages
    for j=i+1:nImages
        ind = ind + 1;
        allvijs = cell2mat(vijs_new_cell{i,j});
        vijs = reshape(allvijs,9,[]);
        vijs = vijs';
        IDX = kmeans(vijs,2);
        class1 = find(IDX==1);
        vijs_averaged(:,:,ind) = mean(allvijs(:,:,class1),3);
    end
end

sign_ij_J = cryo_sync3n_Jsync_power_method(vijs_averaged,1,0);

vijs_out = zeros(size(vijs_averaged));

for i=1:numel(sign_ij_J)
    if sign_ij_J(i) == 1
        vijs_out(:,:,i) = vijs_averaged(:,:,i);
    else
        vijs_out(:,:,i) = J*vijs_averaged(:,:,i)*J;
    end
end

sign_J_dist = histc(sign_ij_J,[-1,1])/numel(sign_ij_J);
log_message('global J-sync dist=[%.2f %.2f]',sign_J_dist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Syncronizing viis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At this stage all vij where i!=j are synchronized at this stage. Thus,
% since for any i and for any j it holds that 
% $vi*vi^{T}vi*vj^{T} = vi*vj^{T}$, and similarly 
% $Jvi*vi^{T}JJvi*vj^{T}J = Jvi*vj^{T}J$ we can independently check for any
% estimate vii whether it it j-conjugated or not.

viis_out = zeros(size(viis));
sign_ii_J = zeros(1,nImages);
for i=1:nImages
    err = 0;
    err_J = 0;
    Rii = viis(:,:,i);
    Rii_J = J*Rii*J;
    
    for j=1:i-1
        ind = uppertri_ijtoind(j,i,nImages);
        Rji = vijs_out(:,:,ind);
        
        err = err + norm(Rji*Rii-Rji,'fro');
        err_J = err_J + norm(Rji*Rii_J-Rji,'fro');
    end
    
    for j=i+1:nImages
        ind = uppertri_ijtoind(i,j,nImages);
        Rij = vijs_out(:,:,ind);
        
        err = err + norm(Rii*Rij-Rij,'fro');
        err_J = err_J + norm(Rii_J*Rij-Rij,'fro');
    end
    
    if err < err_J
        sign_ii_J(i) = 1;
    else
        sign_ii_J(i) = -1;
    end
end

for i=1:nImages
    Rii = viis(:,:,i);
    
    if sign_ii_J(i) == 1
        viis_out(:,:,i) = Rii;
    else %i.e., sign_J_diag(i) == -1
        viis_out(:,:,i) = J*Rii*J;
    end
end

sign_J_diag_dist = histc(sign_ii_J,[-1,1])/numel(sign_ii_J);
log_message('Global J-synch diag dist=[%.2f %.2f]',sign_J_diag_dist);


end
%
% global J;
% % sign_J encodes the partition.
% V_out = zeros(size(V));
% for i=1:nImages
%     for j=i+1:nImages
%         ind = uppertri_ijtoind(i,j,nImages);
%         vij = V((i-1)*3+1:i*3,(j-1)*3+1:j*3);
%         if sign_J(ind) == 1
%             V_out((i-1)*3+1:i*3,(j-1)*3+1:j*3) = vij;
%         else % i.e., sign_J(ind) == -1
%             V_out((i-1)*3+1:i*3,(j-1)*3+1:j*3) = J*vij*J;
%         end
%     end
% end
%
% sign_J_dist = histc(sign_J,[-1,1])/numel(sign_J);
% log_message('outer_sync_dist=[%.2f %.2f]',sign_J_dist);
%
% %synchronize the diagonal blocks. At this stage, the off-diagonal blocks of
% %V_out are either all without a spurios J in them, or all with a sspurious
% %J in them.
% sign_J_diag = zeros(1,nImages);
% for i=1:nImages
%     err = 0;
%     err_J = 0;
%     vii = V((i-1)*3+1:i*3,(i-1)*3+1:i*3);
%     vii_J = J*vii*J;
%
%     for j=1:i-1
%         vji = V_out((j-1)*3+1:j*3,(i-1)*3+1:i*3);
%
%         err = err + norm(vji*vii-vji,'fro');
%         err_J = err_J + norm(vji*vii_J-vji,'fro');
%     end
%
%     for j=i+1:nImages
%         vij = V_out((i-1)*3+1:i*3,(j-1)*3+1:j*3);
%
%         err = err + norm(vii*vij-vij,'fro');
%         err_J = err_J + norm(vii_J*vij-vij,'fro');
%     end
%
%     if err < err_J
%         sign_J_diag(i) = 1;
%     else
%         sign_J_diag(i) = -1;
%     end
% end
%
% for i=1:nImages
%     vii = V((i-1)*3+1:i*3,(i-1)*3+1:i*3);
%     if sign_J_diag(i) == 1
%         V_out((i-1)*3+1:i*3,(i-1)*3+1:i*3) = vii;
%     else %i.e., sign_J_diag(i) == -1
%         V_out((i-1)*3+1:i*3,(i-1)*3+1:i*3) = J*vii*J;
%     end
% end
%
% sign_J_diag_dist = histc(sign_J_diag,[-1,1])/numel(sign_J_diag);
% log_message('outer_sync_diag_dist=[%.2f %.2f]',sign_J_diag_dist);
%
% end