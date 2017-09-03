function [rots,in_plane_rots] = estimate_inplane_rotations4(npf,vis,inplane_rot_res,max_shift,shift_step)
%
% General description
% 
% Input parameters:
%   npf              A 3D array where each image npf(:,:,i) corresponds to the Fourier
%                    transform of projection i.
%   vis              description
%   inplane_rot_res  (Optional) The resolution in degrees of in-plane rotation angle.
%                    E.g., resolution=1 means we look for angles 0,1,2,3 ...
%                    E.g.,  resolution=0.25 means we look for angles 0,0.25,0.5,...
%                    Defualt=1
%   max_shift 
%   shift_step          
% Output parameters:
%   rots       A 3D array of size 3x3xnImages where the i-th slice
%              rots(:,:,3) is equal to Ri (the i-th rotation matrix)

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = 15;
end

if ~exist('inplane_rot_res','var')
    inplane_rot_res = 1;
end

log_message('estimating in-plane rotation angles')

[n_r,n_theta,nImages] = size(npf);

H = zeros(nImages,nImages);

% create all in-plane angular differences using 'inplane_rot_res' as the
% discretization
theta_ij = (0:inplane_rot_res:(360-inplane_rot_res))*pi/180;
n_theta_ij = numel(theta_ij);

% TODO: this should be enforced in the interface
assert(mod(n_theta_ij,5)==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: construct all in-plane rotation matrices R(theta_ij)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cos_theta_ij = cos(theta_ij);
sin_theta_ij = sin(theta_ij);
zrs = zeros(1,n_theta_ij);
ons = ones(1,n_theta_ij);
R_theta_ij = [cos_theta_ij  ; sin_theta_ij; zrs; ...
             -sin_theta_ij  ; cos_theta_ij; zrs; ...
              zrs;            zrs;          ons];

R_theta_ij = reshape(R_theta_ij,3,3,n_theta_ij);

max_corrs = zeros(1,nchoosek(nImages,2));
max_idx_corrs = zeros(1,nchoosek(nImages,2));
counter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: construct all tilde_Ri. These are rotation matrices 
% whose third row is equal to that of Ri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ris_tilde = zeros(3,3,nImages);
for i = 1:nImages
    v_i = vis(:,i).';
    Ris_tilde(:,:,i) = complete_3rdRow_to_rot(v_i);
end

% ignoring dc term
npf_tmp = npf;
npf_tmp(1,:,:) = 0;

npf_normalized = bsxfun(@rdivide,npf_tmp,...
                                 sqrt(sum((abs(npf_tmp)).^2)));

clear npf_tmp;

%precompile the shift phases
shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
[~,nshifts] = size(shift_phases);
npf_i_shifted = zeros(n_r,n_theta,nshifts);
msg = [];
for i = 1:nImages
    
    %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(i,5) == 0
        bs = char(repmat(8,1,numel(msg)));
        fprintf('%s',bs);
        msg = sprintf('k=%3d/%3d',i,nImages);
        fprintf('%s',msg);
    end
    %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    npf_i = npf(:,:,i);
    % get all possible shifted copies of the image
    for s=1:nshifts
        npf_i_shifted(:,:,s) = bsxfun(@times,npf_i,shift_phases(:,s));
    end
    
    npf_i_shifted(1,:) = 0; % ignoring dc term
    % normalize each ray to be of norm 1
    npf_i_shifted = bsxfun(@rdivide,npf_i_shifted,...
        sqrt(sum((abs(npf_i_shifted)).^2)));
    
    Npf_i_shifted = gpuArray(single(npf_i_shifted));
    
    % for performance purposes we precalculate this once per each image i.
    % If it were for readibility it would be better to put this in the inner
    % loop.
    Ri_tilde = Ris_tilde(:,:,i);
    % calculate R_i_tilde*R_theta_ij for all possible R_theta_ijs
    tmp = multiprod(Ri_tilde.',R_theta_ij);
    
    for j = i+1:nImages
        
%         t1 = clock;
        counter = counter+1;
        
        npf_j = npf_normalized(:,:,j);
        
        
        Rj_tilde = Ris_tilde(:,:,j);
        
        % to each of the previous computed matrices we multiply it by Rj_tilde
        Us  = multiprod(tmp,Rj_tilde);
        
        Us = reshape(Us,9,n_theta_ij);
        % extract a common-line index for each possible theta_ij
        c_1 = [-Us(8,:) ;  Us(7,:)]; %[-U(2,3)  U(1,3)];
        c_2 = [ Us(6,:) ; -Us(3,:)]; %[ U(3,2) -U(3,1)];
        
        cijs = clAngles2Ind(c_1,n_theta);
        cjis = clAngles2Ind(c_2,n_theta);
        
        Npf_j = gpuArray(single(npf_j));
        
        % cross-correaltion - so we want to conjugate
        co = bsxfun(@times,Npf_i_shifted(:,cijs,:),conj(Npf_j(:,cjis))); 
        Corrs = sum(co);
        corrs = gather(Corrs);
        
        corrs = reshape(corrs,n_theta_ij/5,5,nshifts);
        
        if nshifts > 1
            % the 1d shifts depend not only on the 2d shifts of each and every
            % image, but also on their common-lines. That is, for a given image, the 1d
            % shift might be different for different lines !
            % In particular, this means that each one of the 4 lines we are considering in a given image
            % may correspond to a DIFFERENT 1d shift. Therefore each line has the dof, so to speak, of choosing its own shift.
            % All in all, this means that we max out the shift BEFORE taking
            % the mean score over all 4 lines.
            corrs = max(real(corrs),[],3);
        end
        % now take the mean score over all 4 lines, and find the maximum
        % quadraple of lines
        [max_corr,max_idx_corr] = max(mean(real(corrs),2));
        max_corrs(counter) = max_corr; % this is only for stats
        max_idx_corrs(counter) = max_idx_corr; % this is only for stats
        theta_diff = inplane_rot_res*(max_idx_corr-1)*pi/180;
        
        H(i,j) = cos(5*theta_diff) + sqrt(-1)*sin(5*theta_diff);
%            
%         %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         t2 = clock;
%         t = etime(t2,t1);
%         bs = char(repmat(8,1,numel(msg)));
%         fprintf('%s',bs);
%         msg = sprintf('k=%3d/%3d k=%3d/%3d t=%7.5f',i,nImages,j,nImages,t);
%         fprintf('%s',msg);
%         
%         %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

fprintf('\n');
% note entry (i,j) corresponds to exp^(i(-theta_i+theta_i)). Therefore to
% construct the hermitian matrix, entry (j,i) is the **conjugate** of entry (i,j)
H = H + H' + eye(nImages); % put 1 on diagonal since : exp^(i*0) = 1

[v, d] = eigs(H, 10, 'lm');
[evals, ind] = sort(diag(d), 'descend');
evect1 = v(:,ind(1)); %zi = exp(-i*2*ai), Ri = Rot(ai)*Ri0

log_message('first 5 eigenvalues=[%.2f %.2f %.2f %.2f %.2f]',evals(1),evals(2),evals(3),evals(4),evals(5));

% figure,
% hist(evals, 40);
% set(gca, 'FontSize', 20);
% title('SO(2) Sync');

% H is rank-1 whose eigenvector is given by (exp^{-2pi theta_1},...,exp^{-2pi theta_n})
in_plane_rots = zeros(3,3,nImages);
for i  = 1:nImages
    zi  = evect1(i);
    zi  = zi/abs(zi);      % rescale so it lies on unit circle
    c   = real(zi^(1/5));  % we have four time the desired angle
    s   = -imag(zi^(1/5)); % we have four time the desired angle
    in_plane_rots(:,:,i) = [c -s 0; s c 0; 0 0 1];
end

log_message('forming rotation matrices');
rots = zeros(3,3,nImages);
for i=1:nImages
    rots(:,:,i) = in_plane_rots(:,:,i)*Ris_tilde(:,:,i);
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