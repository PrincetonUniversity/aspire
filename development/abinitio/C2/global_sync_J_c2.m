function [Rijs_out,Rijgs_out,sign_ij_J] = global_sync_J_c2(Rijs,Rijgs,nImages)
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


log_message('\nGlobal J synchronization');
% partition all off-diagonal blocks into two classes: Those with a spurious
% J in them and those without a spurious J in them.
sz = nchoosek(nImages,2);
assert(size(Rijs,3) == sz);

vijs = zeros(3,3,sz);
for i=1:sz
    Rij  = Rijs(:,:,i);
    Rijg = Rijgs(:,:,i);
    vijs(:,:,i) = (Rij+Rijg)/2;
end

sign_ij_J = cryo_sync3n_Jsync_power_method(vijs,1,0,1);

assert(numel(sign_ij_J) == sz);

Rijs_out  = zeros(size(Rijs));
Rijgs_out = zeros(size(Rijgs));


J = diag([1 1 -1]); % Reflection matrix
for i=1:numel(sign_ij_J)
    if sign_ij_J(i) == 1
        Rijs_out(:,:,i)  = Rijs(:,:,i);
        Rijgs_out(:,:,i) = Rijgs(:,:,i);
    else
        Rijs_out(:,:,i)  = J*Rijs(:,:,i)*J;
        Rijgs_out(:,:,i) = J*Rijgs(:,:,i)*J;
    end
end

sign_J_dist = histc(sign_ij_J,[-1,1])/numel(sign_ij_J);
log_message('global J-sync dist=[%.2f %.2f]',sign_J_dist);

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