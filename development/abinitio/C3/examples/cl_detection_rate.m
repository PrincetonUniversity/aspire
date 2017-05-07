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

clmatrix_diff = bsxfun(@minus,clmatrix_gt,clmatrix);
clmatrix_diff_angle = clmatrix_diff*2*pi./n_theta;
% take absolute cosine because of handedness
% there might be +180 independendt diff for each image which at this stage
% hasn't been taken care yet.
nCorrect = 0;
hand_idx = zeros(1,nchoosek(nImages,2));
for i=1:nImages
    for j=i+1:nImages
        ind = uppertri_ijtoind(i,j,nImages);
        diffs_cij = clmatrix_diff_angle(i,j,:);
        diffs_cji = clmatrix_diff_angle(j,i,:);
        min_diff1 = min(acos(cos(diffs_cij))    + acos(cos(diffs_cji)));
        min_diff2 = min(acos(cos(diffs_cij+pi)) + acos(cos(diffs_cji+pi)));
        if min_diff1 < min_diff2
            min_diff = min_diff1;
            hand_idx(ind) = 1;
        else
            min_diff = min_diff2;
            hand_idx(ind) = 2;
        end
        if min_diff < 2*angle_tol_err
            nCorrect  = nCorrect+1;
            clmatrix_correct(i,j) = 1;
            clmatrix_correct(j,i) = 1;
        end
    end
end

cl_dist = histc(hand_idx,1:2)/numel(hand_idx);
detec_rate = nCorrect/(nImages*(nImages-1)/2);
log_message('common lines detection rate=%.2f%%',detec_rate*100);
log_message('cl_J_dist=[%.2f %.2f]',cl_dist);

end

function clmatrix_gt = find_cl_gt(n_theta,refq)


nImages = size(refq,2);
clmatrix_gt = zeros(nImages,nImages,4);

g = [cosd(120) -sind(120) 0; ...
     sind(120)  cosd(120) 0; ...
     0                 0  1]; % rotation matrix of 120 degress around z-axis

gs = zeros(3,3,3);
for s=0:2
    gs(:,:,s+1) = g^s;
end

for i=1:nImages
    for j=i+1:nImages
        Ri = q_to_rot(refq(:,i))';
        Rj = q_to_rot(refq(:,j))';
        for s=0:2
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