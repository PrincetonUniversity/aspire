function [vijs_out,sign_ij_J] = global_sync_J_Dn(vijs, nImages)
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
%                  the corresponding outer-product +-vi*vj^{T} between the
%                  third rows of matrices Ri and Rj. Each such estimate
%                  might have a spurious J independently of other estimates
%   nImages        The number of images (or rotations) to estimate (N).
% Output parameters:
%   vijs_out     A 3x3xn_choose_2 array where each slice holds an estimate for
%                the corresponding outer-product vi*vj^{T} 
%   sign_ij_J    An array of length n-choose-2 where each entry equals 1
%                or -1 depending whether the corresponding input 
%                estimate vivj had a spurious J or not.
%
% Based on global_sync_J (without vii synchronization).
% Modified by Elad Eatah, May 2021
%


log_message('Global J synchronization');
% partition all off-diagonal blocks into two classes: Those with a spurious
% J in them and those without a spurious J in them.
assert(size(vijs,3) == nchoosek(nImages,2));

sign_ij_J = cryo_syncDn_Jsync_power_method(vijs,1,0);

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
end