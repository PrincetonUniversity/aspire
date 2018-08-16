function [im_avg_est,im_avg_est_orig,log_lik,opt_latent,outlier_ims_inds] = em_class_avg_updated(images,init_avg_image,nIters,ang_jump,max_shift,shift_jump,nScaleTicks,remove_outliers,verbose,gt_data)
% Expectation-Maximization (EM) algorithm to find the image that best fits the input noisy images, 
% assuming each image has was rotated, shifted, and scaled.
% grid to the PSWF expansion coefficients.
%   Input:  images:          A 3D array of noisy images. The third dimension corresponds to the image index.
%           init_avg_image:  (Optional)  The average image to initialize
%                            the em algorithm. Supposedly, all other images are a rotated,
%                            shifted, and scaled version thereof. Default:
%                            the first image in 'images'.
%           nIters:          (Optional) The number of iterations to run the em algorithm. Defualt: 10
%           ang_jump:        (Optional) the resolution of rotation angles for the em to asses. That is, the em will consider the angles
%                             1:ang_jump:360. Default:1
%           max_shift:       (Optional) the maximum number of shifts in pixels for the em to assess in both x and y directions. Default:5. 
%                             provided basis functions correspond to the set of parameters.
%          shift_jump:       (Optional) the resolution of shifts for the em
%                             to asses. That is, the em will consider the shifts
%                             -max_shift:shift_jump:max_shift. Default:1.
%                             Note: it must be that max_shift is a multiple of shift_jump so that shift=0 is considered
%          nScaleTicks       (Optional) the number of different scales for
%                            the em algorithm to use. Defualt:10
%          remove_outliers   (Optional) if true, then a second pass of the
%                             em algorithm will be performed on images that
%                             are not ouliers. Default =false
%          verbose:         (Optional) 0:no output messages, 1:some
%                           intermidiate messages, 2:detailed messages and
%                           debug images. Default=0
%   Output: 
%
%       im_avg_est:          The average image (maximum likelhood over all images over all latent variables)
%       log_lik:               : a (nImages,nIters) array holding the log likelhood of each image for each iteration
%       opt_latent:           : a struct of containing the peak value for
%                               each latent variable for each image

if ~exist('init_avg_image','var') || isempty(init_avg_image)
    init_avg_image = images(:,:,1);
    log_message('initializing initial average image to be the first image');
end

if ~exist('nIters','var')
    nIters = 10;
end

if ~exist('ang_jump','var')
    ang_jump = 1;
end

if ~exist('max_shift','var')
    max_shift = 5;
end

if ~exist('shift_jump','var')
    shift_jump = 1;
end

if ~exist('nScaleTicks','var')
    nScaleTicks = 10;
end

if ~exist('remove_outliers','var')
    remove_outliers = false;
else
    outliers_precent_removal = 10;
end

if ~exist('verbose','var')
    verbose = 0;
end

if ~exist('gt_data','var') || isempty(gt_data)
    gt_data = [];
end

if mod(max_shift,shift_jump) ~= 0
    error('max_shift=%d is not an integral multiple of shift_jump=%d',max_shift,shift_jump);
end

if ~isempty(gt_data)
    assert(isfield(gt_data,'thetas_gt') && isfield(gt_data,'shifts_gt') && isfield(gt_data,'scales_gt'));
end

if verbose >= 1; log_message('#iterations=%d,angualr-jump=%d,max shift=%d,shift-jump=%d,#scales=%d',nIters,ang_jump,max_shift,shift_jump,nScaleTicks); end

%% general parameters
if size(images,1) ~= size(images,2)
    error('The input images must be squares');
end

if size(init_avg_image,1) ~= size(init_avg_image,2)
    error('The initial average image must be square');
end

if size(init_avg_image,1) ~= size(images,1)
    error('The initial average image must be the same size as images');
end

opt_im_size = 65; % the algorithm support images up to size of 129 but it is slow
if size(images,1) > opt_im_size
    is_downsample = true;
else
    is_downsample = false;
end

if is_downsample
    if verbose >= 1; log_message('downsampling the input image to be of size 65 for the em algorithm to run fast'); end
    images_orig         = images;
    init_avg_image_orig = init_avg_image;
    images         = cryo_downsample(images_orig,65,1);
    init_avg_image = cryo_downsample(init_avg_image_orig,65);
end

im_size = size(images,1);
nImages = size(images,3);

%% Prolates-related parameters
L = floor(im_size/2); % Image resolution: (2L+1)x(2L+1) samples on a Cartesian grid
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values stand for greater oversampling
% c = beta*pi*L;  % Bandl imit 
T = 1e+1;       % Truncation parameter

%% mask backgrounds so that we may then assume that sigma=1 for all images
init_avg_image = cryo_mask(init_avg_image);
% [init_avg_image,init_img_mean_bg,init_img_sd_bg] = my_cryo_normalize_background(init_avg_image);
[images,mean_bg_ims,sd_bg_ims] = cryo_normalize_background(images,floor(im_size/2),0);

%% estimate the snr to get a possible range of scales
snr_est   = cryo_estimate_snr(images);
if snr_est <= 0
    snr_est = 10^-4;
end
est_scale = sqrt(snr_est*mean(sd_bg_ims)^2);

%% em-related parameters
em_params.scales = linspace(0.8*est_scale,1.2*est_scale,nScaleTicks);
em_params.thetas = 1:ang_jump:360;

% em_params.shifts = -max_shift:max_shift;
em_params.max_shift = max_shift;
em_params.shift_jump = shift_jump;
em_params.shifts = -max_shift:shift_jump:max_shift;

%% express the images using prolates (extract the expansion coefficients)
if verbose >= 1; log_message('Obtaining all image coefficients while evaluating PSWFs'); end
[c_ims, PSWF_Nn_p] = my_pswf_t_f(images, L, beta, T, []);  % Obtain coefficients while evaluating PSWFs
[c_avg, ~]         = my_pswf_t_f(init_avg_image, L, beta, T, PSWF_Nn_p); 

c_avg      = gpuArray(c_avg);
nRots      = numel(em_params.thetas);
numel_Nn   = size(PSWF_Nn_p.samples,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    EM Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% precompute the rotation phases
if verbose >= 2; log_message('precomputing phases'); end
phases = zeros(numel(PSWF_Nn_p.unq_ang_freqs),nRots);
for t=1:nRots 
    theta = em_params.thetas(t);
    phases(:,t) =  exp(-sqrt(-1)*PSWF_Nn_p.unq_ang_freqs.*theta*2*pi/360);
end

%% precompute the expansion coefficients of each image for each possible rotation
c_ims_rot = zeros([numel_Nn,nRots,nImages]);
phases_full = phases(PSWF_Nn_p.Iunq_ang_freqs,:);
for i=1:nImages
    c_ims_rot(:,:,i) = bsxfun(@times,c_ims(:,i),phases_full);
end
clear phases_full;
% the coefficients that correspond to the negative frequencies are 
% conjugate equal to those of the positive frequencies, so we shall upload
% to the gpu only those
c_ims_rot = c_ims_rot(PSWF_Nn_p.non_neg_freqs,:,:); 
c_ims_rot = gpuArray(c_ims_rot);

%% precompute the non-cross term in the e-step that does not change.
if verbose >= 2; log_message('precomputing constant terms'); end
const_terms = precompute_const_terms(c_ims,mean_bg_ims,sd_bg_ims,im_size,L,beta,T,PSWF_Nn_p);
const_terms = structfun(@gpuArray,const_terms,'UniformOutput',false);
%%
log_lik = cell(1,2);

%%
% These arrays are only needed for stats (simulated data where ground-truth is available)
simulation_data.scales_err  = [];
simulation_data.shift_x_err = [];
simulation_data.shift_y_err = [];

outlier_ims_inds = [];
im_avg_est_prev = init_avg_image;
for round = 1:2
    log_lik{round} = zeros(nImages,nIters); % holds the likilhood of the data given the current parameters 
    for iter=1:nIters
        
        if verbose >=1; log_message('iter %d',iter); end
        %% E-step
        if verbose >=1; log_message('E-step'); end
        [Posteriors,log_lik_per_image] = e_step(c_avg,c_ims_rot,PSWF_Nn_p,const_terms,em_params,mean_bg_ims,sd_bg_ims,verbose);
        log_lik{round}(:,iter) = log_lik_per_image;
        if verbose >=1; log_message('iter %d: log likelihhod=%.2f',iter,sum(log_lik_per_image)); end
        
        %% M-step
        if verbose >=1; log_message('M-step'); end
        c_avg = m_step(Posteriors,c_ims,PSWF_Nn_p,const_terms,em_params,sd_bg_ims,phases,verbose);
        
        %% debug code
        if verbose >= 2
            %% debug code when ground truth data is available
            if ~isempty(gt_data)
                simulation_data = plot_latent_vars_stats(Posteriors,em_params,gt_data,simulation_data);
            end
            im_avg_est = my_pswf_t_b(gather(c_avg),PSWF_Nn_p);
            plot_images(init_avg_image,im_avg_est_prev,im_avg_est);
            im_avg_est_prev = im_avg_est;
        end
    end
    if round == 1 && remove_outliers %maximum two rounds 
        [~,II] = sort(log_lik{round}(:,end),'ascend');
        outlier_ims_inds = II(1:floor(outliers_precent_removal/100*size(images,3)));
        
        Posteriors(:,:,:,outlier_ims_inds) = [];
        images(:,:,outlier_ims_inds) = [];
        c_ims_rot(:,:,outlier_ims_inds) = [];
        
        c_ims(:,outlier_ims_inds) = [];
        mean_bg_ims(outlier_ims_inds) = [];
        sd_bg_ims(outlier_ims_inds) = [];
        const_terms = precompute_const_terms(c_ims,mean_bg_ims,sd_bg_ims,im_size,L,beta,T,PSWF_Nn_p);
        const_terms = structfun(@gpuArray,const_terms,'UniformOutput',false);
        
        if is_downsample
            images_orig(:,:,outlier_ims_inds) = [];
        end
        
        if ~isempty(gt_data)
            gt_data.thetas_gt(outlier_ims_inds)   = [];
            gt_data.shifts_gt(:,outlier_ims_inds) = [];
            gt_data.scales_gt(outlier_ims_inds)   = [];
        end
        nImages = size(images,3);
    else
        break;
    end
end
im_avg_est = real(my_pswf_t_b(gather(c_avg),PSWF_Nn_p));

%% in case the original image was too large, the em algorithm operated on the downsampled imaged. 
%% So we make a final m-step pass on the original images using the latest posteriors 
if is_downsample
    clear PSWF_Nn_p; % too heavy to leave on the gpu
    im_avg_est_orig = do_one_pass_orig_images(Posteriors,images_orig,images,em_params,beta,T,verbose); 
    plot_images(init_avg_image_orig,im_avg_est_orig,im_avg_est_orig);
else
    im_avg_est_orig = im_avg_est;
end

%% find the mode of each latent variable for each image
opt_latent = compute_opt_latent_vals(Posteriors,em_params);

end

function im_avg_est_orig = do_one_pass_orig_images(Posteriors,images_orig,images,em_params,beta,T,verbose)

if verbose >= 1; log_message('applying a single m-step on original images using latest posteriors'); end
im_size_small = size(images,1);
im_size_orig = size(images_orig,1);
L_orig = floor(im_size_orig/2); % Image resolution: (2L+1)x(2L+1) samples on a Cartesian grid
[images_orig,mean_bg_ims_orig,sd_bg_ims_orig] = cryo_normalize_background(images_orig);

%% shift each image according to the mode of the shift using the latest posterier array

%step 1: marginilize over latent variables which are not shifts, and
% find the optimal shift per image
if verbose >= 1; log_message('finding optimal shift per image using latest posteriors'); end
nShifts_1d = numel(em_params.shifts);
[nScales,nRots,nShifts_2d,nImages] = size(Posteriors);
[~,opt_shift_ind_per_image] = max(sum(reshape(Posteriors,[nScales*nRots,nShifts_2d,nImages])));

Posteriors_opt_shift = zeros([nScales,nRots,1,nImages],'gpuArray'); % a single shift. Need a 4D array to be able to call the m-step as usual
opt_shift_ind_per_image = squeeze(opt_shift_ind_per_image);
for i=1:nImages
    post_i = Posteriors(:,:,opt_shift_ind_per_image(i),i);
    post_i = post_i/sum(post_i(:));
    Posteriors_opt_shift(:,:,1,i) = post_i;
end


if verbose >= 1; log_message('shifting each image by the optimal shift found'); end
[yys,xxs] = ind2sub([nShifts_1d,nShifts_1d],opt_shift_ind_per_image);

opt_shifts_x = em_params.shifts(xxs);
opt_shifts_y = em_params.shifts(yys);

size_ratio = round(im_size_orig/im_size_small); % TODO: do interpolation
opt_shifts_x = size_ratio*opt_shifts_x;
opt_shifts_y = size_ratio*opt_shifts_y;
for i=1:nImages
    images_orig(:,:,i) = circshift(images_orig(:,:,i),[opt_shifts_x(i),opt_shifts_y(i)]);
end

% step 2: Obtain coefficients while evaluating PSWFs
if verbose >= 1; log_message('obtaining coefficients of shifted images'); end
[c_ims_orig, PSWF_Nn_p_orig] = my_pswf_t_f(images_orig, L_orig, beta, T, []);

%     PSWF_Nn_p_orig.Psis = gather(PSWF_Nn_p_orig.Psis); % We might not have enough gpu memeory for the psis
if verbose >= 1; log_message('precomputing constant terms'); end
const_terms_orig = precompute_const_terms(c_ims_orig,mean_bg_ims_orig,sd_bg_ims_orig,im_size_orig,L_orig,beta,T,PSWF_Nn_p_orig);

if verbose >= 2; log_message('precomputing phases'); end
phases_orig = zeros(numel(PSWF_Nn_p_orig.unq_ang_freqs),nRots);
for t=1:nRots
    theta = em_params.thetas(t);
    phases_orig(:,t) =  exp(-sqrt(-1)*PSWF_Nn_p_orig.unq_ang_freqs.*theta*2*pi/360);
end

% since we shifted each image using its mode we need not consider
% shifts any more
em_params.max_shift = 0;
em_params.shifts = 0;
c_avg_orig = m_step(Posteriors_opt_shift,c_ims_orig,PSWF_Nn_p_orig,const_terms_orig,em_params,sd_bg_ims_orig,phases_orig,verbose);

im_avg_est_orig = real(my_pswf_t_b(gather(c_avg_orig),PSWF_Nn_p_orig));

end

function const_terms = precompute_const_terms(c_ims,mean_bg_ims,sd_bg_ims,im_sz,L_orig,beta,T,PSWF_Nn_p_orig)

nImages = size(c_ims,2);
% we need the all ones image in order to accomodate for the additive term due to normalization
[c_all_ones_im, ~] = my_pswf_t_f(ones(im_sz,im_sz), L_orig, beta, T, PSWF_Nn_p_orig);
const_terms.c_all_ones_im = gpuArray(c_all_ones_im);

const_terms.anni = zeros(1,nImages);
for i=1:nImages
    const_terms.anni(i) = norm(c_ims(:,i))^2;
end

const_terms.cnn = gpuArray(norm(const_terms.c_all_ones_im)^2*((mean_bg_ims./sd_bg_ims).^2));

const_terms.c_additive_term = bsxfun(@times,const_terms.c_all_ones_im,mean_bg_ims./sd_bg_ims);

end

function opt_latent = compute_opt_latent_vals(Posteriors,em_params)

[~,~,~,nImages] = size(Posteriors);

nShifts_1d = numel(em_params.shifts);

opt_latent.rots     = zeros(1,nImages);
opt_latent.shifts_x = zeros(1,nImages);
opt_latent.shifts_y = zeros(1,nImages);
opt_latent.scales   = zeros(1,nImages);
for i=1:nImages
    om_i = Posteriors(:,:,:,i);
    
    [~,opt_scale_ind] = max(sum(sum(om_i,2),3));
    [~,opt_rot_ind]   = max(sum(sum(om_i,1),3));
    [~,opt_shift_ind] = max(sum(sum(om_i,1),2));
    
    opt_latent.scales(i) = em_params.scales(opt_scale_ind);
    opt_latent.rots(i) = em_params.thetas(opt_rot_ind);
    
    [yy,xx] = ind2sub([nShifts_1d,nShifts_1d],opt_shift_ind);
    opt_latent.shifts_x(i) = em_params.shifts(xx);
    opt_latent.shifts_y(i) = em_params.shifts(yy);
    
end

end

function A_shift = calc_A_shift(Psis,shift_x,shift_y,PSWF_Nn_p)

numel_Nn = size(Psis,3);
neg_feq_inds  = PSWF_Nn_p.neg_feq_inds;
zero_feq_inds = PSWF_Nn_p.zero_feq_inds;
pos_feq_inds  = PSWF_Nn_p.pos_feq_inds;

non_neg_freqs = PSWF_Nn_p.non_neg_freqs;

Psis_shftd = circshift(Psis(:,:,non_neg_freqs),[shift_x,shift_y]);
Psis_shftd = bsxfun(@times,PSWF_Nn_p.points_inside_the_circle,Psis_shftd);

Psis       = reshape(Psis,      [],numel_Nn);
Psis_shftd = reshape(Psis_shftd,[],numel(non_neg_freqs));

A_shift(:,non_neg_freqs) = Psis'*Psis_shftd; % we need the conjugation by design
A_shift(zero_feq_inds,neg_feq_inds) = conj(A_shift(zero_feq_inds,pos_feq_inds));
A_shift(pos_feq_inds,neg_feq_inds)  = conj(A_shift(neg_feq_inds,pos_feq_inds));
A_shift(neg_feq_inds,neg_feq_inds)  = conj(A_shift(pos_feq_inds,pos_feq_inds));

end


function points_inside_the_circle = get_points_inside_the_circle(L)

x_1d_grid = -L:1:L;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= L);

end



function [coeffs, PSWF_Nn_p] = my_pswf_t_f(images, L, beta, T, PSWF_Nn_p)

realFlag = 0;
if isempty(PSWF_Nn_p)
    [coeffs, PSWF_Nn_p] = pswf_t_f(images, L, beta, T, realFlag, PSWF_Nn_p);
    
    PSWF_Nn_p.samples = [ conj(PSWF_Nn_p.samples(:,PSWF_Nn_p.ang_freq>0)) PSWF_Nn_p.samples];    
    PSWF_Nn_p.ang_freq = [-PSWF_Nn_p.ang_freq(PSWF_Nn_p.ang_freq>0) ; PSWF_Nn_p.ang_freq];
    
    PSWF_Nn_p.neg_feq_inds = find(PSWF_Nn_p.ang_freq<0);
    PSWF_Nn_p.zero_feq_inds = find(PSWF_Nn_p.ang_freq==0);
    PSWF_Nn_p.pos_feq_inds = find(PSWF_Nn_p.ang_freq>0);
    PSWF_Nn_p.non_neg_freqs = [PSWF_Nn_p.zero_feq_inds ; PSWF_Nn_p.pos_feq_inds]';
    PSWF_Nn_p.points_inside_the_circle = get_points_inside_the_circle(PSWF_Nn_p.L);
    [unq_ang_freqs,~,Iunq_ang_freqs] = unique(PSWF_Nn_p.ang_freq); %PSWF_Nn_p.ang_freq = unq_ang_freqs(Iunq_ang_freqs)
    PSWF_Nn_p.unq_ang_freqs = unq_ang_freqs;
    PSWF_Nn_p.Iunq_ang_freqs = Iunq_ang_freqs;
    
    numel_Nn = size(PSWF_Nn_p.samples,2);
    L = PSWF_Nn_p.L;
    
    if L > floor(129/2)
        Psis = zeros([2*L+1,2*L+1,numel_Nn]);
    else
        Psis = zeros([2*L+1,2*L+1,numel_Nn],'gpuArray');
    end
    
    for nn=1:numel_Nn
        PSWF_sample = PSWF_Nn_p.samples(:,nn);
        Psi_nn = zeros(2*L+1,2*L+1);
        Psi_nn(PSWF_Nn_p.points_inside_the_circle) = PSWF_sample;
        Psis(:,:,nn) = Psi_nn;
    end
    PSWF_Nn_p.Psis = Psis;
    
else
    
    PSWF_Nn_p.samples(:,PSWF_Nn_p.neg_feq_inds) = [];
    PSWF_Nn_p.ang_freq(PSWF_Nn_p.neg_feq_inds) = [];
    
    [coeffs, PSWF_Nn_p] = pswf_t_f(images, L, beta, T, realFlag, PSWF_Nn_p);
    
    PSWF_Nn_p.samples = [ conj(PSWF_Nn_p.samples(:,PSWF_Nn_p.ang_freq>0)) PSWF_Nn_p.samples];    
    PSWF_Nn_p.ang_freq = [-PSWF_Nn_p.ang_freq(PSWF_Nn_p.ang_freq>0) ; PSWF_Nn_p.ang_freq];    
end

end

function images = my_pswf_t_b( coeffs, PSWF_Nn_p)

realFlag = 0;
PSWF_Nn_p.samples(:,PSWF_Nn_p.neg_feq_inds) = [];
PSWF_Nn_p.ang_freq(PSWF_Nn_p.neg_feq_inds) = [];

images = pswf_t_b( coeffs, PSWF_Nn_p, realFlag );

PSWF_Nn_p.samples = [ conj(PSWF_Nn_p.samples(:,PSWF_Nn_p.ang_freq>0)) PSWF_Nn_p.samples];
PSWF_Nn_p.ang_freq = [-PSWF_Nn_p.ang_freq(PSWF_Nn_p.ang_freq>0) ; PSWF_Nn_p.ang_freq];    

end


function [Posteriors,log_lik_per_image] = e_step(c_avg,c_ims_rot,PSWF_Nn_p,const_terms,em_params,mean_bg_ims,sd_bg_ims,verbose)

nShifts_1d = numel(em_params.shifts);
nShifts_2d = nShifts_1d^2; % both x and y directions
nScales    = numel(em_params.scales);
nRots      = numel(em_params.thetas);
numel_Nn   = size(PSWF_Nn_p.samples,2);
nImages = numel(mean_bg_ims);
non_neg_freqs = PSWF_Nn_p.non_neg_freqs;
anni = const_terms.anni;
cnn  = const_terms.cnn;
c_all_ones_im = const_terms.c_all_ones_im;

% The posterior array
Posteriors = zeros(nScales,nRots,nShifts_2d,nImages,'gpuArray');

%% compute the terms that do not depend on the shifts
ann_const      =  gpuArray(((em_params.scales.'*(1./sd_bg_ims)).^2)*norm(c_avg)^2);
cross_cnn_ann  =  gpuArray(em_params.scales.'*(mean_bg_ims./(sd_bg_ims.^2))*2*real(c_avg'*c_all_ones_im));
ann_const_cross_cnn_anns = ann_const + cross_cnn_ann;

const_elms = bsxfun(@plus,ann_const_cross_cnn_anns,anni + cnn);

c_ims_i_rot_full = zeros(numel_Nn,nRots,'gpuArray'); % will hold the exapnasion coeffs for each possible rotation of a given image
if verbose >= 2; printProgressBarHeader; end
counter = 0;
% The heaviest operation is to calculate the the 'A' matrix of a given
% shift. But since A(shift) = A(-shift)' we need only iterate over half
% of all shifts and manully compute the others by taking the conjugate
% transpose
for shift_x=em_params.shifts
    for shift_y=em_params.shifts
        
        if shift_y < shift_x; continue; end
        
        counter  = counter + 1; % just for debugging
        
        A_shift = calc_A_shift(PSWF_Nn_p.Psis,shift_x,shift_y,PSWF_Nn_p);
        
        %% pre calculate stuff that are image invariant and depend only on the shift
        
        % stuff resulated to (shift_x,shift_y)
        tmp1_shift = c_all_ones_im' * A_shift;
        tmp2_shift = c_avg' * A_shift;
        
        % stuff resulated to (-shift_x,-shift_y)
        A_inv_shift = A_shift';
        tmp1_inv_shift = c_all_ones_im' * A_inv_shift;
        tmp2_inv_shift = c_avg' * A_inv_shift;
        
        %translate from 2-d shift to 1-d shift index
        subinds_I = ([ shift_y -shift_y] + em_params.max_shift)/em_params.shift_jump + 1;
        subinds_J = ([ shift_x -shift_x] + em_params.max_shift)/em_params.shift_jump + 1;
        inds = sub2ind([nShifts_1d,nShifts_1d],subinds_I,subinds_J);
        
        for i=1:nImages
            
            if verbose >= 2; progressTic((counter-1)*nImages+i,nImages*(nShifts_1d + 0.5*(nShifts_1d)^2)); end
            
            % in the gpu we only kept the coeffs of the non-negative
            % frequencies. Now we need to take all freqs
            c_ims_i_rot_full(non_neg_freqs,:) = c_ims_rot(:,:,i);
            c_ims_i_rot_full(PSWF_Nn_p.neg_feq_inds,:) = conj(c_ims_i_rot_full(PSWF_Nn_p.pos_feq_inds,:));
            
            % calculate the two cross terms
            cross_anni_cnn = mean_bg_ims(i)/sd_bg_ims(i) * 2*real(tmp1_shift * c_ims_i_rot_full);
            cross_anni_ann = em_params.scales.'/sd_bg_ims(i) * 2*real(tmp2_shift * c_ims_i_rot_full);
            
            % write down the log likelihood
            Posteriors(:,:,inds(1),i) = cross_anni_ann - bsxfun(@plus,const_elms(:,i),cross_anni_cnn);
            
            if shift_x ~= shift_y % then we also popolute the the inverse shift to the omega
                cross_anni_cnn_minus = mean_bg_ims(i)/sd_bg_ims(i) * 2*real(tmp1_inv_shift * c_ims_i_rot_full);
                cross_anni_ann_minus = em_params.scales.'/sd_bg_ims(i) * 2*real(tmp2_inv_shift * c_ims_i_rot_full);
                
                Posteriors(:,:,inds(2),i) = cross_anni_ann_minus - bsxfun(@plus,const_elms(:,i),cross_anni_cnn_minus);
            end      
        end
    end
end

%%
% Change from log space to regular space and compute 
% a. the liklihod per image, 
% b. the likelhood of all images and, 
% c. the posterior
log_lik_per_image = zeros(1,nImages); 
for i=1:nImages
    
    omega_i   = Posteriors(:,:,:,i);
    max_omega = max(omega_i(:));
    
    %we multiply by the largest number to avoid numerical problems (implicetely this also multiplies the denominator)
    omega_i = exp(omega_i - max_omega);
    
    log_lik_per_image(i) = gather(max_omega + log(sum(omega_i(:)))); % for stats
    
    % normalize the log likelhood to obtain the posteriors
    Posteriors(:,:,:,i) = omega_i/sum(omega_i(:));
end

end


function c_avg = m_step(Posteriors,c_ims,PSWF_Nn_p,const_terms,em_params,sd_bg_ims,phases,verbose)


[nScales,nRots,~,nImages] = size(Posteriors);
numel_Nn = size(PSWF_Nn_p.Psis,3);
nShifts_1d = numel(em_params.shifts);
non_neg_freqs = PSWF_Nn_p.non_neg_freqs;
n_unq_ang_freqs = numel(PSWF_Nn_p.unq_ang_freqs);
Iunq_ang_freqs = PSWF_Nn_p.Iunq_ang_freqs;

% calculate the normalization number
c = bsxfun(@times,bsxfun(@times,Posteriors,reshape(em_params.scales,[nScales,1,1,1])),reshape(1./sd_bg_ims,[1,1,1,nImages]));
c = sum(c(:));
Phases = gpuArray(phases.'); % note the transpose. Maybe fix later on
C_ims         = gpuArray(c_ims);
c_avg         = zeros(numel_Nn,1,'gpuArray');
W_shifts_marg = zeros(n_unq_ang_freqs,nImages,'gpuArray');
if verbose >= 2; printProgressBarHeader; end
counter = 0;
for shift_x=em_params.shifts
    for shift_y=em_params.shifts
        
        if shift_y < shift_x; continue; end
        
        %% move from 2-d shift index to 1d shift index
        counter  = counter + 1;
        subinds_I = ([ shift_y -shift_y] + em_params.max_shift)/em_params.shift_jump + 1;
        subinds_J = ([ shift_x -shift_x] + em_params.max_shift)/em_params.shift_jump + 1;
        
        inds = sub2ind([nShifts_1d,nShifts_1d],subinds_I,subinds_J);
        
        if verbose >= 2; progressTic(counter,nShifts_1d + 0.5*(nShifts_1d)^2); end
        
        if shift_x == 0 && shift_y == 0
            A_shift = eye(numel_Nn,numel_Nn);
        else
            A_shift = calc_A_shift(PSWF_Nn_p.Psis,shift_x,shift_y,PSWF_Nn_p); 
        end
        
        A_inv_shift = A_shift';
        
        % we **know** that the coeffs of the negative freqs are conjugate equal
        % to the coeffs of the positive freqs, so no point in calculating explicetely both
        A_shift       = A_shift(non_neg_freqs,:);
        A_inv_shift   = A_inv_shift(non_neg_freqs,:);
        
        % marginilize out the numbers that correspond to scales and rotations
        W = zeros(n_unq_ang_freqs,nImages,'gpuArray');
        for i=1:nImages
            W(:,i) = sum(Posteriors(:,:,inds(1),i)*Phases);
%             W_tmp = bsxfun(@times,reshape(Phases,[1,nRots,n_unq_ang_freqs]),Posteriors(:,:,inds(1),i));
%             W(:,i) = sum(reshape(W_tmp,nScales*nRots,n_unq_ang_freqs));
        end
        % update all non-negs freqs at once
        c_avg(non_neg_freqs) = c_avg(non_neg_freqs) + sum(A_shift*(W(Iunq_ang_freqs,:).*C_ims),2);
        
        % the additive term is shift invariant as far as the posteriors
        % are concerned. In particular shifts should be marginilized as
        % well
        W_shifts_marg = W_shifts_marg + W;
        
        if shift_x ~= shift_y % then we update also the inverse shift
            W_minus = zeros(n_unq_ang_freqs,nImages,'gpuArray');
            for i=1:nImages
                W_minus(:,i) = sum(Posteriors(:,:,inds(2),i)*Phases);
%                 W_minus_tmp = bsxfun(@times,reshape(Phases,[1,nRots,n_unq_ang_freqs]),Posteriors(:,:,inds(2),i));
%                 W_minus(:,i) = sum(reshape(W_minus_tmp,nScales*nRots,n_unq_ang_freqs));
            end
            
            c_avg(non_neg_freqs) = c_avg(non_neg_freqs) + sum(A_inv_shift*(W_minus(Iunq_ang_freqs,:).*C_ims),2);
            W_shifts_marg = W_shifts_marg + W_minus;
        end
    end
end

%% update the coeffs using with respect to the additive term
A_no_shift = eye(numel_Nn,numel_Nn); %calc_A_shift(PSWF_Nn_p.Psis,0,0,PSWF_Nn_p);
A_no_shift = A_no_shift(non_neg_freqs,:);
c_avg(non_neg_freqs) = c_avg(non_neg_freqs) + sum(A_no_shift*(W_shifts_marg(Iunq_ang_freqs,:).*const_terms.c_additive_term),2);
%% populate the coeffs of the negative frequencies
c_avg(PSWF_Nn_p.neg_feq_inds) = conj(c_avg(PSWF_Nn_p.pos_feq_inds));
%% normalize
c_avg = c_avg/c;

clear Exp_minus_intheta;
clear C_ims;

end


function simulation_data = plot_latent_vars_stats(Posteriors,em_params,gt_data,simulation_data)

prev_scales_err  = simulation_data.scales_err;
prev_shift_x_err = simulation_data.shift_x_err;
prev_shift_y_err = simulation_data.shift_y_err;

[nScales,nRots,nShifts_2d,nImages] = size(Posteriors);
nShifts_1d = numel(em_params.shifts);

opt_rots     = zeros(1,nImages);
opt_shifts_x = zeros(1,nImages);
opt_shifts_y = zeros(1,nImages);
opt_scales   = zeros(1,nImages);

figure(2);
for k=1:nImages
    om_k = Posteriors(:,:,:,k);
    [~,ii] = max(om_k(:));
    [opt_scale,opt_rot,opt_shift] = ind2sub([nScales,nRots,nShifts_2d],ii);
    
    opt_scales(k) = em_params.scales(opt_scale);
    
    opt_rots(k) = em_params.thetas(opt_rot);
    
    [yy,xx] = ind2sub([nShifts_1d,nShifts_1d],opt_shift);
    
    opt_shifts_x(k) = em_params.shifts(xx);
    opt_shifts_y(k) = em_params.shifts(yy);
end

%
%         err(iter) = norm(im_avg_est(:)-im_avg_gt(:))/norm(im_avg_gt(:));
%         log_message('iter %d: error=%.4f',iter,err(iter));

subplot(2,2,1); plot(acosd(abs(cosd(opt_rots-gt_data.thetas_gt))),'.');    axis([1 nImages 0 360]);    title('Opt rots err diff');

scales_err = opt_scales-gt_data.scales_gt;
% to generate each image we shift the average image. But in the
% algorithm we apply the inverse shift to the image so that it
% matched the average image. Hence the err difference we need to
% apply minus
shift_x_err = -opt_shifts_x-gt_data.shifts_gt(1,:);
shift_y_err = -opt_shifts_y-gt_data.shifts_gt(2,:);

if ~isempty(prev_scales_err) && ~isempty(prev_shift_x_err) && ~isempty(prev_shift_y_err)
    subplot(2,2,2); plot(scales_err,'.'); hold on;  plot(prev_scales_err,'r.'); hold off; title('Opt scale err diff');
    subplot(2,2,3); plot(shift_x_err,'.'); hold on;  plot(prev_shift_x_err,'r.'); hold off; axis([1 nImages -2*em_params.shifts(end)-1 2*em_params.shifts(end)+1]); title('Opt shift-x err diff');
    subplot(2,2,4); plot(shift_y_err,'.'); hold on;  plot(prev_shift_y_err,'r.'); hold off; axis([1 nImages -2*em_params.shifts(end)-1 2*em_params.shifts(end)+1]); title('Opt shift-y err diff');
else
    subplot(2,2,2); plot(scales_err,'.'); title('Opt scale err diff');
    subplot(2,2,3); plot(shift_x_err,'.'); axis([1 nImages -2*em_params.shifts(end)-1 2*em_params.shifts(end)+1]); title('Opt shift-x err diff');
    subplot(2,2,4); plot(shift_y_err,'.'); axis([1 nImages -2*em_params.shifts(end)-1 2*em_params.shifts(end)+1]); title('Opt shift-y err diff');
end

drawnow;

simulation_data.scales_err  = scales_err;
simulation_data.shift_x_err = shift_x_err;
simulation_data.shift_y_err = shift_y_err;


end

function plot_images(init_avg_image,im_avg_est_prev,im_avg_est)

h = figure(1);

subplot(1,4,1); imagesc(real(init_avg_image));  colormap(gray); axis off; title('First');
subplot(1,4,2); imagesc(real(im_avg_est_prev)); colormap(gray); axis off; title('Prev');
subplot(1,4,3); imagesc(real(im_avg_est));      colormap(gray); axis off; title('Current');
subplot(1,4,4); imagesc(real(im_avg_est-im_avg_est_prev)); colormap(gray); axis off; title('Diff');

truesize(h);

drawnow;

end