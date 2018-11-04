function [detec_rate,clmatrix_correct] = cl_detection_rate_c2(clmatrix,n_theta,min_dist_cls,refq)
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


angle_tol_err = 5/180*pi; % how much angular deviation we allow for a common-line to have
nImages = size(clmatrix,1);
clmatrix_correct = zeros(size(clmatrix));
% clmatrix_gt is a n*n*2 matrix representing the four pairs of common-lines between each two images
clmatrix_gt = find_cl_gt(n_theta,refq);
foo(clmatrix_gt);
% calculate the diff of FIRST common-line in each image against the two
% lines of ground-truth
clmatrix_diff_11 = (clmatrix_gt(:,:,1) - clmatrix(:,:,1))*2*pi./n_theta;
clmatrix_diff_12 = (clmatrix_gt(:,:,1) - clmatrix(:,:,2))*2*pi./n_theta;
clmatrix_diff_21 = (clmatrix_gt(:,:,2) - clmatrix(:,:,1))*2*pi./n_theta;
clmatrix_diff_22 = (clmatrix_gt(:,:,2) - clmatrix(:,:,2))*2*pi./n_theta;

% take absolute cosine because of handedness
% there might be +180 independendt diff for each image which at this stage
% hasn't been taken care yet.
% hand_idx = zeros(nchoosek(nImages,2),2);
cls_flip = zeros(1,nchoosek(nImages,2));
hands    = zeros(2,nchoosek(nImages,2));
inds_one_found = [];
clFound = 0;
nzeroFound = 0;
nOneFound = 0;
ntwoFound = 0;
for i=1:nImages
    for j=i+1:nImages
        ind = uppertri_ijtoind(i,j,nImages);
        
        diffs_cij_11 = acos(cos(clmatrix_diff_11(i,j)));
        diffs_cij_12 = acos(cos(clmatrix_diff_12(i,j)));
        diffs_cij_21 = acos(cos(clmatrix_diff_21(i,j)));
        diffs_cij_22 = acos(cos(clmatrix_diff_22(i,j)));
        
        diffs_cji_11 = acos(cos(clmatrix_diff_11(j,i)));
        diffs_cji_12 = acos(cos(clmatrix_diff_12(j,i)));
        diffs_cji_21 = acos(cos(clmatrix_diff_21(j,i)));
        diffs_cji_22 = acos(cos(clmatrix_diff_22(j,i)));
        
        
        diffs_cij_11_J = acos(cos(clmatrix_diff_11(i,j)+pi));
        diffs_cij_12_J = acos(cos(clmatrix_diff_12(i,j)+pi));
        diffs_cij_21_J = acos(cos(clmatrix_diff_21(i,j)+pi));
        diffs_cij_22_J = acos(cos(clmatrix_diff_22(i,j)+pi));
        
        diffs_cji_11_J = acos(cos(clmatrix_diff_11(j,i)+pi));
        diffs_cji_12_J = acos(cos(clmatrix_diff_12(j,i)+pi));
        diffs_cji_21_J = acos(cos(clmatrix_diff_21(j,i)+pi));
        diffs_cji_22_J = acos(cos(clmatrix_diff_22(j,i)+pi));
                      
        is_cls_found = zeros(1,4);
        hand_found   = zeros(1,4);
        if diffs_cij_11 < angle_tol_err && diffs_cji_11 < angle_tol_err
            clmatrix_correct(i,j,1) = clmatrix(i,j,1);
            clmatrix_correct(j,i,1) = clmatrix(j,i,1);
            is_cls_found(1) = 1;
            hand_found(1) = 0;
        elseif diffs_cij_11_J < angle_tol_err && diffs_cji_11_J < angle_tol_err
            clmatrix_correct(i,j,1) = clmatrix(i,j,1);
            clmatrix_correct(j,i,1) = clmatrix(j,i,1);
            is_cls_found(1) = 1;
            hand_found(1) = 1;
        end
        
        if diffs_cij_22 < angle_tol_err && diffs_cji_22 < angle_tol_err
            clmatrix_correct(i,j,2) = clmatrix(i,j,2);
            clmatrix_correct(j,i,2) = clmatrix(j,i,2);
            is_cls_found(2) = 1;
            hand_found(2) = 0;
        elseif diffs_cij_22_J < angle_tol_err && diffs_cji_22_J < angle_tol_err
            clmatrix_correct(i,j,2) = clmatrix(i,j,2);
            clmatrix_correct(j,i,2) = clmatrix(j,i,2);
            is_cls_found(2) = 1;
            hand_found(2) = 1;
        end
        
        if diffs_cij_12 < angle_tol_err && diffs_cji_12 < angle_tol_err
            clmatrix_correct(i,j,2) = clmatrix(i,j,2);
            clmatrix_correct(j,i,2) = clmatrix(j,i,2);
            is_cls_found(3) = 1;
            hand_found(3) = 0;
        elseif diffs_cij_12_J < angle_tol_err && diffs_cji_12_J < angle_tol_err
            clmatrix_correct(i,j,2) = clmatrix(i,j,2);
            clmatrix_correct(j,i,2) = clmatrix(j,i,2);
            is_cls_found(3) = 1;
            hand_found(3) = 1;
        end
        
        if diffs_cij_21 < angle_tol_err && diffs_cji_21 < angle_tol_err
            clmatrix_correct(i,j,1) = clmatrix(i,j,1);
            clmatrix_correct(j,i,1) = clmatrix(j,i,1);
            is_cls_found(4) = 1;
            hand_found(4) = 0;
        elseif diffs_cij_21_J < angle_tol_err && diffs_cji_21_J < angle_tol_err
            clmatrix_correct(i,j,1) = clmatrix(i,j,1);
            clmatrix_correct(j,i,1) = clmatrix(j,i,1);
            is_cls_found(4) = 1;
            hand_found(4) = 1;
        end
                
        [nFound,ii] = max([is_cls_found(1)+is_cls_found(2), is_cls_found(3)+is_cls_found(4)]); 
        clFound = clFound + nFound;
        cls_flip(ind) = ii;
        if ii==1
            hands(:,ind) = [hand_found(1); hand_found(2)];
        else
            hands(:,ind) = [hand_found(3); hand_found(4)];
        end
        if nFound == 0
            nzeroFound = nzeroFound +1;
        elseif nFound == 1
            nOneFound = nOneFound +1;
            inds_one_found(:,end+1) = [i;j];
        else
            ntwoFound = ntwoFound + 1;
        end
    end
end

per = analyze_one_found(inds_one_found,clmatrix_gt,min_dist_cls,n_theta);

detec_rate = clFound/(2*nchoosek(nImages,2));
log_message('common lines detection rate=%.2f%%',detec_rate*100);

detec_rate_no_cl = nzeroFound/(nchoosek(nImages,2));
log_message('common lines detection rate zero cl=%.2f%%',detec_rate_no_cl*100);

detec_rate_single_cl = nOneFound/(nchoosek(nImages,2));
log_message('common lines detection rate one cl=%.2f%%',detec_rate_single_cl*100);

detec_rate_two_cl = ntwoFound/(nchoosek(nImages,2));
log_message('common lines detection rate two cl=%.2f%%',detec_rate_two_cl*100);

cl_J_dist = histc(hands(:),0:1)/numel(hands(:));
log_message('cl_J_dist=[%.2f %.2f]',cl_J_dist);

cl_flip_dist = histc(cls_flip(:),1:2)/numel(cls_flip);
log_message('cl_flip_dist=[%.2f %.2f]',cl_flip_dist)

end

function clmatrix_gt = find_cl_gt(n_theta,refq)


nImages = size(refq,2);
clmatrix_gt = zeros(nImages,nImages,2);

g = diag([-1 -1 1]);


for i=1:nImages
    for j=i+1:nImages
        Ri = q_to_rot(refq(:,i))';
        Rj = q_to_rot(refq(:,j))';
        
        % first common-lines
        U = Ri.'*Rj;
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
        clmatrix_gt(i,j,1) = idx1;
        clmatrix_gt(j,i,1) = idx2;
        
        
        % second common-lines
        U = Ri.'*g*Rj;
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
        clmatrix_gt(i,j,2) = idx1;
        clmatrix_gt(j,i,2) = idx2;     
    end
end


end


function per = analyze_one_found(inds_one_found,clmatrix_gt,min_dist_cls,n_theta)

min_dist_cls = min_dist_cls*pi/180;

nPairs = size(inds_one_found,2);

nclose_pairs = 0;
for k = 1:nPairs
    i = inds_one_found(1,k);
    j = inds_one_found(2,k);
    
    diff_ij = acos(cos((clmatrix_gt(i,j,1) - clmatrix_gt(i,j,2))*2*pi./n_theta));
    diff_ji = acos(cos((clmatrix_gt(j,i,1) - clmatrix_gt(j,i,2))*2*pi./n_theta));
    
    if diff_ij < min_dist_cls && diff_ji < min_dist_cls
        nclose_pairs = nclose_pairs + 1; 
    end
end

per = 100*nclose_pairs/nPairs;

end


function foo(clmatrix_gt)

nImages = size(clmatrix_gt,1);
nPairs = nchoosek(nImages,2);
diffs = 10:10:180;
stats = zeros(numel(diffs),nPairs);

for i=1:nImages
    for j=i+1:nImages
        
        cij_1 = clmatrix_gt(i,j,1);
        cij_2 = clmatrix_gt(i,j,2);
        cji_1 = clmatrix_gt(j,i,1);
        cji_2 = clmatrix_gt(j,i,2);
        
        
        ind = uppertri_ijtoind(i,j,nImages);
        for k=1:numel(diffs)
            diff = diffs(k);
            stats(k,ind) = ((acosd(cosd(cij_1 - cij_2)) < diff) && (acosd(cosd(cji_1 - cji_2))) < diff);
%             stats(k,ind) = ((acosd(cosd(cij_1 - cij_2)) < diff) || (acosd(cosd(cji_1 - cji_2))) < diff);
        end
    end
end

percent_per_diff = squeeze(100*sum(stats,2)/nPairs);
% figure; 
% plot(percent_per_diff(:,1)); hold on; 
% plot(percent_per_diff(:,2));
end