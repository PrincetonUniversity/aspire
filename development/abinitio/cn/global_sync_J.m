function [vijs_out,viis_out,sign_ij_J,sign_ii_J] = global_sync_J(vijs,viis)
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
nImages = size(viis,3);
assert(size(vijs,3) == nchoosek(nImages,2));

sign_ij_J = cryo_sync3n_Jsync_power_method(vijs,1,0);

assert(numel(sign_ij_J) == nchoosek(nImages,2));

vijs_out = zeros(size(vijs));

J = diag([1 1 -1]); % Reflection matrix
for i=1:numel(sign_ij_J)
    if sign_ij_J(i) == 1
        vijs_out(:,:,i) = vijs(:,:,i);
    else
        vijs_out(:,:,i) = J*vijs(:,:,i)*J;
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