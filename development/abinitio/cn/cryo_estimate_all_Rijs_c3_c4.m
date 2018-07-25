function Rijs = cryo_estimate_all_Rijs_c3_c4(n_symm,clmatrix,n_theta,refq)
%
% Estimate a single relative rotation Rij between images i and j. For n_symm==3 (i.e., c3) the
% estimate may correspond to either RiRj, RigRj, or Rig^{2}Rj.
% For n_symm==4 (i.e., c4) the estimate may correspond to either RiRj, RigRj, or Rig^{2}Rj, or Rig^{3}Rj
% 
% Input parameters:
%   n_symm      Either 3 (for c_3) or 4 (for c_4)
%   clmatrix   An n-by-n matrix (where n represents the number of images).
%              The (i,j)-th entry is the index of one of the common lines
%              in image i with image j
%   n_theta    The angular discretization of lines in each image. E.g.,
%              n_theta=360 means that there are 360 lines in each image
%   refq       (Optional) A 4-by-n table. The i-th column represent the quaternion of
%              that corresponds to the rotation matrix of the i-th image
%
% Output parameters:
%   Rijs       A 3x3xn_choose_2 array holding the estimates for the
%              relative orientations


Rijs = cryo_sync3n_estimate_all_Rijs(clmatrix, n_theta);

%%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('refq','var') && ~isempty(refq)
    
g = [cosd(360/n_symm) -sind(360/n_symm) 0; ...
     sind(360/n_symm)  cosd(360/n_symm) 0; ...
     0                 0  1]; % rotation matrix of 120 or 90 degress around z-axis
    
    J = diag([1 1 -1]); % Reflection matrix
    
    nImages = size(refq,2);    
    errs = zeros(1,nchoosek(nImages,2));
    %precompile g,g^2
    gs = zeros(3,3,n_symm);
    for s=0:n_symm-1
        gs(:,:,s+1) = g^s;
    end
    
    for i=1:nImages
        for j=i+1:nImages
            
            ind = uppertri_ijtoind(i,j,nImages);
            Rij = Rijs(:,:,ind);
            
            Ri_gt = q_to_rot(refq(:,i))';
            Rj_gt = q_to_rot(refq(:,j))';
            
            errs_tmp = zeros(1,3);
            for s=1:n_symm
                Rij_gt = Ri_gt.'*gs(:,:,s)*Rj_gt;
                % we are oblivious to a possible J conjugation at this moment
                errs_tmp(s) = min([norm(Rij-Rij_gt,'fro'),norm(J*Rij*J-Rij_gt,'fro')]);
            end
            % we are oblivious to which relative rotation we actualy have
            % in hand, so take the optimal
            errs(ind) = min(errs_tmp);
        end
    end
    
    mse = mean(errs.^2);
    log_message('MSE of Rij''s=%.2f',mse);
end

%%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

end