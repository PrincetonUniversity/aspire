function Rijs = cryo_c2_estimate_all_Rijs(clmatrix,n_theta,refq)
%
% Estimate a single relative rotation Rij between images i and j. The
% estimate may correspond to either RiRj, RigRj, Rig^{2}Rj, or Rig^{3}Rj
% 
% Input parameters:
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
    
    g = [-1  0 0; ...
          0 -1 0; ...
          0  0 1]; % rotation matrix of 180 degress around z-axis
    
    J = diag([1 1 -1]); % Reflection matrix
    
    nImages = size(refq,2);    
    errs = zeros(1,nchoosek(nImages,2));
    
    for i=1:nImages
        for j=i+1:nImages
            
            ind = uppertri_ijtoind(i,j,nImages);
            Rij = Rijs(:,:,ind);
            
            Ri_gt = q_to_rot(refq(:,i))';
            Rj_gt = q_to_rot(refq(:,j))';
            
            Rij_gt  = Ri_gt.'*Rj_gt;
            Rijg_gt = Ri_gt.'*g*Rj_gt;
            % we are oblivious to a possible J conjugation at this moment
            err1 = min([norm(Rij-Rij_gt,'fro'),norm(J*Rij*J-Rij_gt,'fro')]);
            err2 = min([norm(Rij-Rijg_gt,'fro'),norm(J*Rij*J-Rijg_gt,'fro')]);
            
            % we are oblivious to which relative rotation we actualy have
            % in hand, so take the optimal
            errs(ind) = min(err1,err2);
        end
    end
    
    mse = mean(errs.^2);
    log_message('MSE of Rij''s=%.2f',mse);
end

%%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

end