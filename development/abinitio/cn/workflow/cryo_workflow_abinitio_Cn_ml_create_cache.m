function cache_full_file_name  = cryo_workflow_abinitio_Cn_ml_create_cache

message = 'Load existing viewing directions inds cache file? ';
is_load_cache_file = multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if is_load_cache_file
    cache_full_file_name = '';
    while isempty(cache_full_file_name)
        cache_full_file_name = fmtinput('Enter full path viewing directions inds cache file: ','','%s');
        if exist(cache_full_file_name,'file')~=2
            fprintf('viewing directions inds cache file %s does not exist\n', cache_full_file_name);
            cache_full_file_name = '';
        end
    end
    is_valid = validate_inds_cache_file(cache_full_file_name);
    if ~is_valid
        fprintf('Viewing directions inds cache is not valid. Please create a new one');
        fprintf('Aborting...\n');
        return;
    end
else
    message = 'Create viewing directions inds cache file?';
    do_create = multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
    if do_create == 1
        message = 'This should take about 10 mins. Contintue?';
        do_continue = multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
        if do_continue
            cache_dir = fmtinput('Enter full path of output directory for cache file: ','','%s');
            if ~exist(cache_dir,'dir') % Do we need to create directory?
                message = 'Output directory does not exist. Create?';
                do_create = multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
                if do_create == 1
                    mkdir(cache_dir);
                else
                    fprintf('Aborting...\n');
                    return;
                end
            else % Check if directory empty
                listing = dir(cache_dir);
                if numel(listing)>2
                    message = 'Output directory not empty. Continue?';
                    do_continue = multichoice_question(message,{'Y','N'},[ 1, 0],'N');
                    if do_continue == 0
                        fprintf('Aborting...\n');
                        return;
                    end
                end
            end
            
            n_Points_sphere = fmtinput('Number of viewing directions candidates? ',1000,'%d');
            n_theta = fmtinput('Number of samples in the angular direction (n_theta) ',360,'%d');
            inplane_rot_res_cache = 1;
            log_message('Using %d in-plane degrees resolution',inplane_rot_res_cache);
                        
            fprintf('Calculating viewing direction candidates. Please wait...\n');
            [Ris_tilde,R_theta_ijs] = generate_cand_rots(n_Points_sphere,inplane_rot_res_cache,false,[]);
            cijs_inds               = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,false,[]);
            
            cache_full_file_name = sprintf('ml_cn_cache_points%d_ntheta%d_res%d.mat',n_Points_sphere,n_theta,inplane_rot_res_cache);
            cache_full_file_name = fullfile(cache_dir,cache_full_file_name);
            log_message('\nSaving inds to cache file=%s',cache_full_file_name);
            save(cache_full_file_name,'cijs_inds','Ris_tilde','R_theta_ijs','n_theta','-v7.3');
            log_message('\nCache file has been successfuly saved!\n');
            is_valid = validate_inds_cache_file(cache_full_file_name);
            if ~is_valid
                fprintf('Somethin went wrong. Viewing directions inds cache is not valid. Try again');
                fprintf('Aborting...\n');
                return;
            end
        else
            fprintf('Aborting...\n');
            return;
        end
    else
        fprintf('Aborting...\n');
        return;
    end
end

end


function [Ris_tilde,R_theta_ijs] = generate_cand_rots(nPoints_sphere,inplane_rot_res,is_use_gt_in_cands,refq)

vis = generate_vis(nPoints_sphere,is_use_gt_in_cands,refq);
nRisTilde = size(vis,1);
Ris_tilde  = zeros(3,3,nRisTilde);
for i=1:nRisTilde
    vi = vis(i,:);
    Ris_tilde(:,:,i) = complete_3rdRow_to_rot(vi);
end

theta_ij = (0:inplane_rot_res:(360-inplane_rot_res))*pi/180;
n_theta_ij = numel(theta_ij);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: construct all in-plane rotation matrices R(theta_ij)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cos_theta_ij = cos(theta_ij);
sin_theta_ij = sin(theta_ij);
zrs = zeros(1,n_theta_ij);
ons = ones(1,n_theta_ij);
R_theta_ijs = [cos_theta_ij  ; sin_theta_ij; zrs; ...
    -sin_theta_ij  ; cos_theta_ij; zrs; ...
    zrs;            zrs;          ons];

R_theta_ijs = reshape(R_theta_ijs,3,3,n_theta_ij);

end


function vis = generate_vis(nPoints_sphere,is_use_gt,refq)

if ~is_use_gt
    vis = randn(nPoints_sphere,3);
    for i=1:nPoints_sphere
        % normalize (TODO: what happens if norm is too small. need to discard)
        vi = vis(i,:);
        vi = vi/norm(vi);
%         while (abs( acosd(vi(3)) - 90 ) < 10)
%             vi = randn(1,3);
%             vi = vi/norm(vi);
%         end
        vis(i,:) = vi;
    end
    
    nImages = size(refq,2);
    vis_gt = zeros(nImages,3);
    for i=1:nImages
        rot_gt = q_to_rot(refq(:,i)).';
        vis_gt(i,:) = rot_gt(3,:);
    end
    
    if false %TODO: extract as parameter
        test_coverage(vis.',vis_gt.','red','blue');
    end
else
    nImages = size(refq,2);
    vis_gt = zeros(nImages,3);
    for i=1:nImages
        rot_gt = q_to_rot(refq(:,i))';
        vis_gt(i,:) = rot_gt(3,:);
        
    end
    
    vis_pert = zeros(nImages,3);
    for i=1:nImages
        % add some pertrubation to each row
        noise = randn(1,3);
        noise = tand(5)*noise./norm(noise);
        vis_pert(i,:) = vis_gt(i,:) + noise;
        vis_pert(i,:) = vis_pert(i,:)./norm(vis_pert(i,:));
    end
    
    vis_bad = zeros(2*nImages,3);
    for i=1:size(vis_bad,1)
        tmp = randn(3,3); [Q,~] = qr(tmp);
        vis_bad(i,:) = Q(3,:);
    end
        
    vis = [vis_gt ; vis_pert; vis_bad];
end

end

function R = complete_3rdRow_to_rot(r3)
%
% Constructs a rotation matrix whose third row is equal to a given row vector
% 
% Input parameters:
%   r3         A 1x3 vector of norm 1
% Output parameters:
%   R          A rotation matrix whose third row is equal to r3

assert(abs(norm(r3)-1)<1e-5); 
% handle the case that the third row coincides with the z-axis
if norm(r3-[0,0,1]) < 1e-5
    R = eye(3);
    return;
end

% tmp is non-zero since r3 does not coincide with the z-axis
tmp = sqrt(r3(1)^2 + r3(2)^2);
% contruct an orthogonal row vector of norm 1
r1 = [r3(2)/tmp -r3(1)/tmp 0];
% construct r2 so that r3=r1xr2
r2 = [r3(1)*r3(3)/tmp r3(2)*r3(3)/tmp -tmp];

R = [r1; r2; r3];

end



function test_coverage(vis,vis_gt,color1,color2)

figure;
plot_view_dirs(vis,color1);
hold on;
plot_view_dirs(vis_gt,color2);
% pause;

end



function [cijs_inds,Ris_tilde,R_theta_ijs,n_theta] = compute_cls_inds(Ris_tilde,R_theta_ijs,n_theta,is_save_inds_to_cache,file_cache_name)

nRisTilde = size(Ris_tilde,3);
n_theta_ij = size(R_theta_ijs,3);

cijs = zeros(nRisTilde,nRisTilde,n_theta_ij,2,'uint16');
msg = [];
for i=1:nRisTilde
    Ri_tilde = Ris_tilde(:,:,i);
    % calculate R_i_tilde*R_theta_ij for all possible R_theta_ijs
    tmp = multiprod(Ri_tilde.',R_theta_ijs);
    for j=1:nRisTilde
        t1 = clock;
        Rj_tilde = Ris_tilde(:,:,j);
        Rij = multiprod(tmp,Rj_tilde);
        
        Rij = reshape(Rij,9,n_theta_ij);
        % extract a common-line index for each possible theta_ij
        c_1 = [-Rij(8,:) ;  Rij(7,:)]; %[-U(2,3)  U(1,3)];
        c_2 = [ Rij(6,:) ; -Rij(3,:)]; %[ U(3,2) -U(3,1)];
        
        c_1s = clAngles2Ind(c_1,n_theta);
        c_2s = clAngles2Ind(c_2,n_theta);
        
        inds_tmp = find(c_1s > n_theta/2);
        c_1s(inds_tmp) = c_1s(inds_tmp) - n_theta/2;
        c_2s(inds_tmp) = c_2s(inds_tmp) + n_theta/2;
        c_2s(inds_tmp) = mod(c_2s(inds_tmp)-1,n_theta)+1;
        
        cijs(i,j,:,1) = c_1s;
        cijs(i,j,:,2) = c_2s;
        
        %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t2 = clock;
        t = etime(t2,t1);
        bs = char(repmat(8,1,numel(msg)));
        fprintf('%s',bs);
        msg = sprintf('i=%3d/%3d j=%3d/%3d t=%7.5f',i,nRisTilde,j,nRisTilde,t);
        fprintf('%s',msg);
        %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

cijs_inds = uint16(sub2ind([n_theta/2,n_theta],cijs(:,:,:,1),cijs(:,:,:,2)));
clear cijs;

if is_save_inds_to_cache
    log_message('saving inds to cache file=%s',file_cache_name);
    save(file_cache_name,'cijs_inds','Ris_tilde','R_theta_ijs','n_theta','-v7.3');
end
end