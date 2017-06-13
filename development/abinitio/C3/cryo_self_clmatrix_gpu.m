function [sclmatrix,correlations,shifts] = cryo_self_clmatrix_gpu(npf,max_shift,shift_step,refq)
% Input parameters:
%   npf                     A 3D array where each image npf(:,:,i) corresponds to the Fourier
%                           transform of projection i.
%   max_shift               The maximum spatial shift that each image is
%                           assumed to have. Default:15
%   shift_step              (Optional) Default:0.5
%   refq                    (Optional) A 2d table where the i-th column is the
%                           quaternion that corresponds to the beaming direction of
%                           the i-th image.
% Output parameters:
%   sclmatrix               A 2xnImages table where the i-th column holds
%                           the indexes of the first and third
%                           self-common-line in image i. 
%   correlations            An array of length nImages whose i-th entry
%                           holds the correlation between the two self common-lines that were found
%                           in image i
%   shifts                  An array of length nImages whose i-th entry
%                           holds the shift found for image i

log_message('detecting self-common-lines');

[n_r,n_theta,nImages] = size(npf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: detect self-common-lines for all images
% Step 2: detect all equator images
% Step 3: detect self-common-lines for all equator images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters for the angle between self-common-lines. In theory it is
% [60,180] but all 180 apart lines are perfectly correlated so we set it to
% be smaller.
min_angle_diff = 60*pi/180;
max_angle_diff = 165*pi/180;

% the self-common-line matrix holds per image two indeces that represent
% the two self common-lines in the image
sclmatrix = zeros(2,nImages);
correlations = zeros(1,nImages);
shifts = zeros(1,nImages);

% the angular difference between each two self common lines in a given
% image is [90,180], so create a mask.
[X,Y] = meshgrid(1:n_theta/2,1:n_theta);
diff = Y-X;
unsigned_angle_diff = acos(cos(diff.*2*pi./n_theta));

good_diffs = unsigned_angle_diff > min_angle_diff & ...
    unsigned_angle_diff < max_angle_diff;

shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
[n_r_,nshifts] = size(shift_phases);
assert(n_r == n_r_);
g_shift_phases = gpuArray(single(shift_phases));

msg = [];
for i=1:nImages
    
    t1 = clock;
    
    npf_i = npf(:,:,i);
    % each image is conjugate symmetric so need only consider half of the lines
    pi_half = npf_i(:,1:n_theta/2);
    
    g_npf_i   = gpuArray(single(npf_i));
    g_pi_half = gpuArray(single(pi_half));
    
    % generate all shifted copies of the image
    g_pi_half_shifted = zeros([size(pi_half),nshifts],'gpuArray');
    for s=1:nshifts
        g_pi_half_shifted(:,:,s) = bsxfun(@times,g_pi_half,g_shift_phases(:,s));
    end
    
    g_pi_half_shifted = reshape(g_pi_half_shifted,n_r,n_theta/2*nshifts);
    
    % ignoring dc-term
    g_npf_i(1,:) = 0;
    g_pi_half_shifted(1,:) = 0;
    
    % nomalize each ray to be norm 1
    norms = sqrt(sum((abs(g_npf_i)).^2));
    g_npf_i = bsxfun(@rdivide,g_npf_i,norms);
    
    % nomalize each ray to be norm 1
    norms = sqrt(sum((abs(g_pi_half_shifted)).^2));
    g_pi_half_shifted = bsxfun(@rdivide,g_pi_half_shifted,norms);
    
    Corr = g_npf_i.'*g_pi_half_shifted;
    corr = gather(Corr);
    corr = reshape(corr,[n_theta,n_theta/2,nshifts]);
    
    corr = bsxfun(@times,corr,good_diffs);
    
    [correlation,idx] = max(real(corr(:)));
    
    [sc1, scl3, shift] = ind2sub([n_theta, n_theta/2, nshifts],idx);
    sclmatrix(:,i) = [sc1, scl3]';
    correlations(i) = correlation;
    shifts(i) = shift; %TODO: need to "translate" shift 2*max_shift+1 = ...
    
    %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2 = clock;
    t = etime(t2,t1);
    bs = char(repmat(8,1,numel(msg)));
    fprintf('%s',bs);
    msg = sprintf('k=%3d/%3d t=%7.5f',i,nImages,t);
    fprintf('%s',msg);
    %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

log_message('\nself-common-lines correlation (median)=%.2f',median(correlations));

% analysis against ground-truth
if exist('refq','var') && ~isempty(refq)
    scl_detection_rate(sclmatrix,n_theta,refq);
end

end


function detec_rate = scl_detection_rate(sclmatrix,n_theta,refq)
%
% Calculates the detection rate of self-common-lines. It is invariant to a
% possible J-ambiguity (independently for each image), as well as invariant
% to an RigRi<->Rig^{3}Ri ambiguity (independently for each image)
% 
% Input parameters:
%   sclmatrix  A 2xnImages table where the i-th column holds
%              the indexes of the first and third
%              self-common-line in image i. 
%   n_theta    The angular resolution of the images. E.g.,n_theta=360 
%              means that there are 360 lines in each images 
%   refq       A 4-by-n table. The i-th column represent the quaternion of
%              that corresponds to the rotation matrix of the i-th image
%
% Output parameters:
%   detec_rate The detection rate in [0,1] of self common-lines against the
%               ground truth
              
angle_tol_err = 5/180*pi;
% Two issues:
% 1. DOF: cannot tell the difference between first\third self-common-line
% 2. Ambiguity: handedness - At this stage the handedness is independent for each image
nImages = size(sclmatrix,2);
% sclmatrix_correct = zeros(size(sclmatrix));
sclmatrix_gt = detectScls_gt(n_theta,refq); % clmatrix_gt is a 2*n matrix

sclmatrix_diff1 = sclmatrix_gt - sclmatrix;
sclmatrix_diff2 = sclmatrix_gt - flipud(sclmatrix); % we cannot (and need not) tell the difference between the first scl and third scl
clmatrix_diff1_angle = sclmatrix_diff1*2*pi./n_theta;
clmatrix_diff2_angle = sclmatrix_diff2*2*pi./n_theta;

% 1. cosine is invariant to 2pi.
% 2. abs is invariant to +-pi diff corresponding to J-ambiguity
clmatrix_diff1_angle_mean = mean(acos(abs(cos(clmatrix_diff1_angle))));
clmatrix_diff2_angle_mean = mean(acos(abs(cos(clmatrix_diff2_angle))));

% we need not tell the difference between the first scl and third scl so
% take the min difference
[min_mean_angle_diff, scl_idx] = min([clmatrix_diff1_angle_mean ; clmatrix_diff2_angle_mean]);
correct_idxs = find(min_mean_angle_diff < angle_tol_err);
detec_rate = numel(correct_idxs)/nImages;
% for debug purposes: just making sure that more or less half of time we
% retrieved the first self-cl and half of the time the third self-cl
scl_dist = histc(scl_idx,1:2)/numel(scl_idx);
log_message('self-common lines detection rate=%.2f%%',detec_rate*100);
log_message('scl-distribution=[%.2f %.2f]',scl_dist);


% find the polar angles of viewing directions of offending rotations
rots_gt = zeros(3,3,nImages);
for k = 1:nImages
    rots_gt(:,:,k) = q_to_rot(refq(:,k))';
end

bad_idxs  = find(min_mean_angle_diff >= angle_tol_err);
if ~isempty(bad_idxs)
    bad_polar_angs = acosd(rots_gt(3,3,bad_idxs));
    nbad = min(numel(bad_polar_angs),5);
    log_message('sample of polar angles failing self-cl [%.2f,%.2f,%.2f,%.2f,%.2f]',...
        bad_polar_angs(1:nbad));
end

end

function sclmatrix_gt = detectScls_gt(n_theta,refq)


nImages = size(refq,2);
% we find the first and third self common-lines (c4 symmetry)
sclmatrix_gt = zeros(2,nImages);
% correlations_selfcl  = zeros(1,nImages);
rots_gt = zeros(3,3,nImages);
for i = 1:nImages
    rots_gt(:,:,i) = q_to_rot(refq(:,i))';
end


g = [cosd(120) -sind(120) 0; ...
     sind(120)  cosd(120) 0; ...
     0                 0  1]; % rotation matrix of 120 degress around z-axis


for i=1:nImages
    Ri = rots_gt(:,:,i);
    
    % first self common-lines
    U1 = Ri.'*g*Ri;
    c1=[-U1(2,3)  U1(1,3)]';
    idx1 = clAngles2Ind(c1,n_theta);
    sclmatrix_gt(1,i) = idx1;
    % third self common-lines
    U2 = Ri.'*g^2*Ri;
    c2=[-U2(2,3)  U2(1,3)]';
    idx2 = clAngles2Ind(c2,n_theta);
    sclmatrix_gt(2,i) = idx2;
    
    %npf_k = nonTop_npf(:,:,i);
    %correlations_selfcl(i) = npf_k(:,selfCL_matrix(1,i)).'*npf_k(:,selfCL_matrix(3,i));
end

% if strcmp(params.SCL,'GT') && params.confuse_scl
%     heads = find(round(rand(1,nImages)));
%     sclmatrix_gt(:,heads) = flipud(sclmatrix_gt(:,heads));
% end

end
