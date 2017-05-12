function [detec_rate,clmatrix_correct] = cl_detection_rate(clmatrix,n_theta,refq)
%
% Checks the detection rate of common-lines between 
% images of a c4 symmetric molecule which is invariant to handedness ambiguity. 
% For each pair of images (i,j) the method checks if the
% single pair of common line found is one of the four pairs of common lines
% that are shared by the images.
% 
% Input parameters:
%   clmatrix    A n-by-n table where n represens the number of images.  
%               The (i,j) entry contains the index of a
%               common-line in image i betwen images i and j. 
%   n_theta     The angular resolution of common-lines. It is the number
%               of lines in any given image. E.g., n_theta=360 means that
%               each image has 360 lines
%   refq        A 4-by-n table. The i-th column represent the quaternion of
%               that corresponds to the rotation matrix of the i-th image
%
%
% Output parameters:
%   detec_rate        The detection rate (in [0,1]) of common-lines.
%                      For example, detec_rate=1 means that for each pair
%                      of images the common-line found is one of the four
%                      common-lines between the corresponding images
%
%   clmatrix_correct  A boolean matrix of size n-by-n. the (i,j)-th entry
%                     is equal 1 if the common line found is one of the 
%                     four pairs of common lines between images i and j,
%                     and is 0 otherwise. 

angle_tol_err = 10/180*pi; % how much angular deviation we allow for a common-line to have
nImages = size(clmatrix,1);
clmatrix_correct = zeros(size(clmatrix));
% clmatrix_gt is a n*n*4 matrix representing the four pairs of common-lines between each two images
clmatrix_gt = find_cl_gt(n_theta,refq); 

% calculate the diff of FIRST common-line in each image against the two
% lines of ground-truth
clmatrix_diff_11 = (clmatrix_gt(:,:,1) - clmatrix(:,:,1))*2*pi./n_theta;
clmatrix_diff_12 = (clmatrix_gt(:,:,1) - clmatrix(:,:,2))*2*pi./n_theta;
clmatrix_diff_21 = (clmatrix_gt(:,:,2) - clmatrix(:,:,1))*2*pi./n_theta;
clmatrix_diff_22 = (clmatrix_gt(:,:,2) - clmatrix(:,:,2))*2*pi./n_theta;

% take absolute cosine because of handedness
% there might be +180 independendt diff for each image which at this stage
% hasn't been taken care yet.
nCorrect = 0;
hand_idx = zeros(nchoosek(nImages,2),2);
nclmatrix_correct = zeros(nImages,nImages);
for i=1:nImages
    for j=i+1:nImages
        ind = uppertri_ijtoind(i,j,nImages);
        
        diffs_cij_11 = clmatrix_diff_11(i,j);
        diffs_cij_12 = clmatrix_diff_12(i,j);
        diffs_cij_21 = clmatrix_diff_21(i,j);
        diffs_cij_22 = clmatrix_diff_22(i,j);
        
        diffs_cji_11 = clmatrix_diff_11(j,i);
        diffs_cji_12 = clmatrix_diff_12(j,i);
        diffs_cji_21 = clmatrix_diff_21(j,i);
        diffs_cji_22 = clmatrix_diff_22(j,i);
        
        
        [val_11,hand_11] = min([acos(cos(diffs_cij_11))    + acos(cos(diffs_cji_11)),...
                                acos(cos(diffs_cij_11+pi)) + acos(cos(diffs_cji_11+pi))]);
                      
        
        [val_22,hand_22] = min([acos(cos(diffs_cij_22))    + acos(cos(diffs_cji_22)),...
                                acos(cos(diffs_cij_22+pi)) + acos(cos(diffs_cji_22+pi))]);
                            
        
        [val_12,hand_12] = min([acos(cos(diffs_cij_12))    + acos(cos(diffs_cji_12)),...
                                acos(cos(diffs_cij_12+pi)) + acos(cos(diffs_cji_12+pi))]);
                            
        
        [val_21,hand_21] = min([acos(cos(diffs_cij_21))    + acos(cos(diffs_cji_21)),...
                                acos(cos(diffs_cij_21+pi)) + acos(cos(diffs_cji_21+pi))]);
                  
        
        nCorrect_ij = 0;
        if val_11+val_22 < val_12+val_21
            if val_11 < 2*angle_tol_err
                nCorrect_ij = nCorrect_ij+1;
            end
            if val_22 < 2*angle_tol_err
                nCorrect_ij = nCorrect_ij+1;
            end
            hand_idx(ind,:) = [hand_11 hand_22];
        else
            if val_12 < 2*angle_tol_err
                nCorrect_ij = nCorrect_ij+1;
            end
            if val_21 < 2*angle_tol_err
                nCorrect_ij = nCorrect_ij+1;
            end
            hand_idx(ind,:) = [hand_12 hand_21];
        end
        nclmatrix_correct(i,j) = nCorrect_ij;
        nclmatrix_correct(j,i) = nCorrect_ij;    
        
        nCorrect = nCorrect + nCorrect_ij;
    end
end

cl_J_dist = histc(hand_idx(:),1:2)/numel(hand_idx);
detec_rate = nCorrect/(2*nchoosek(nImages,2));
log_message('common lines detection rate=%.2f%%',detec_rate*100);
log_message('cl_J_dist=[%.2f %.2f]',cl_J_dist);

end

function clmatrix_gt = find_cl_gt(n_theta,refq)


nImages = size(refq,2);
clmatrix_gt = zeros(nImages,nImages,2);

g = [-1  0 0; ...
      0 -1 0; ...
      0  0 1]; % rotation matrix of 180 degress around z-axis

gs = zeros(3,3,2);
for s=0:1
    gs(:,:,s+1) = g^s;
end

for i=1:nImages
    for j=i+1:nImages
        Ri = q_to_rot(refq(:,i))';
        Rj = q_to_rot(refq(:,j))';
        for s=0:1
            U = Ri.'*gs(:,:,s+1)*Rj;
            c1 = [-U(2,3)  U(1,3)]';
            c2 = [ U(3,2) -U(3,1)]';
            
            idx1 = clAngles2Ind(c1,n_theta);
            idx2 = clAngles2Ind(c2,n_theta);
            
%             if strcmp(params_simul.CL,'GT') && params_simul.confuse_cl_J
%                 if round(rand)==1
%                     % j-conjugating amounts at choosing the antipodal
%                     % common-lines in each image
%                     idx1 = mod(idx1+n_theta/2-1,n_theta)+1;
%                     idx2 = mod(idx2+n_theta/2-1,n_theta)+1;
%                 end
%             end
            clmatrix_gt(i,j,s+1) = idx1;
            clmatrix_gt(j,i,s+1) = idx2;
        end
    end
end


end