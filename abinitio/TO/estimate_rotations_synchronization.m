function [rots, H] = estimate_rotations_synchronization(est_rel_rots, u_G, cache_file_name)

% Estimates the rotations using synchronization method.
%
%   Input:
%       est_rel_rots -      n_images x n_images array, containing the
%                           candidates indices forming the relative rotation
%                           Rij = R(:,:,est_rel_rots(i,j))*R(:,:,est_rel_rots(j,i)).'
%       u_G -               J synchronization of the relative rotations.
%       cache_file_name -   contains 'R', the candidates rotation set.
%
%   Output:
%       rots -              3 x 3 x n_images array, rots(:,:,i) is the
%                           rotation of image i.
%       H -                 3 rows synchronization matrices of size 
%                           3*n_images x 3*n_images.
%
%   Written by Adi Shasha January 2021. 


% load the candidates set 'R'.
load(cache_file_name,'R');

n_images = size(est_rel_rots,1);
J = [1 0 0;0 1 0;0 0 -1];


%%% memory allocation and constants %%%

% H(:,:,m) - the m'th synchronization matrix, m=1,2,3.
H = zeros(3*n_images,3*n_images,3);  

% single entry matrices.
e_kl = zeros(3,3,9);
e_kl(:,:,1) = [ 1  0  0; 0  0  0; 0  0  0];      % kl = 11 
e_kl(:,:,2) = [ 0  1  0; 0  0  0; 0  0  0];      % kl = 12   
e_kl(:,:,3) = [ 0  0  1; 0  0  0; 0  0  0];      % kl = 13
e_kl(:,:,4) = [ 0  0  0; 1  0  0; 0  0  0];      % kl = 21 
e_kl(:,:,5) = [ 0  0  0; 0  1  0; 0  0  0];      % kl = 22   
e_kl(:,:,6) = [ 0  0  0; 0  0  1; 0  0  0];      % kl = 23
e_kl(:,:,7) = [ 0  0  0; 0  0  0; 1  0  0];      % kl = 31
e_kl(:,:,8) = [ 0  0  0; 0  0  0; 0  1  0];      % kl = 32
e_kl(:,:,9) = [ 0  0  0; 0  0  0; 0  0  1];      % kl = 33


%%% constructing H %%%

% block (1,2) %
pair_ind = uppertri_ijtoind_vec(1,2,n_images);
H(1:3,4:6,1) = R(:,:,est_rel_rots(1,2))*e_kl(:,:,1)*R(:,:,est_rel_rots(2,1)).';
H(1:3,4:6,2) = R(:,:,est_rel_rots(1,2))*e_kl(:,:,5)*R(:,:,est_rel_rots(2,1)).';
H(1:3,4:6,3) = R(:,:,est_rel_rots(1,2))*e_kl(:,:,9)*R(:,:,est_rel_rots(2,1)).';
if u_G(pair_ind) < 0
    H(1:3,4:6,:) = multi_Jify(H(1:3,4:6,:));
end

% blocks (1,j) j=3,..,n_images %
for p_j=3:n_images
    max_norm_1 = 0;
    max_norm_2 = 0;
    max_norm_3 = 0;
    pair_ind = uppertri_ijtoind_vec(1,p_j,n_images);
    for kl=[1,5,9]
        % H1j = v(k)_1*v(l)_pj.'
        % H(1,2,t) = v(t)_1*v(?)_2.'
        % H(1,2,t).'*H1j = 0 if t != k
        %                = 1 if t == k
        H1j = R(:,:,est_rel_rots(1,p_j))*e_kl(:,:,kl)*R(:,:,est_rel_rots(p_j,1)).';
        if u_G(pair_ind) < 0
            H1j = J*H1j*J;
        end
        norm_1 = norm(H(1:3,4:6,1).'*H1j);
        if norm_1 > max_norm_1
            H(1:3,3*p_j-2:3*p_j,1) = H1j;
            max_norm_1 = norm_1;
        end
        norm_2 = norm(H(1:3,4:6,2).'*H1j);
        if norm_2 > max_norm_2
            H(1:3,3*p_j-2:3*p_j,2) = H1j;
            max_norm_2 = norm_2;
        end
        norm_3 = norm(H(1:3,4:6,3).'*H1j);
        if norm_3 > max_norm_3
            H(1:3,3*p_j-2:3*p_j,3) = H1j;
            max_norm_3 = norm_3;
        end
    end
end

% blocks (i,j) i=2,..,n_images-1, j=3,..,n_images, i<j %
for p_i=2:(n_images-1)
    for p_j=(p_i+1):n_images
        min_norm_1 = 10;
        min_norm_2 = 10;
        min_norm_3 = 10;
        % (v(t)_1*v(k)_pi.').'*(v(t)_1*v(l)_pj.') = 
        % (v(k)_pi*v(t)_1.')*(v(t)_1*v(l)_pj.') = v(k)_pi*v(l)_pj.'
        % where t in {1,2,3} and kl in {1,..,9}
        Hij_1 = H(1:3,3*p_i-2:3*p_i,1).'*H(1:3,3*p_j-2:3*p_j,1);
        Hij_2 = H(1:3,3*p_i-2:3*p_i,2).'*H(1:3,3*p_j-2:3*p_j,2);
        Hij_3 = H(1:3,3*p_i-2:3*p_i,3).'*H(1:3,3*p_j-2:3*p_j,3);
        pair_ind = uppertri_ijtoind_vec(p_i,p_j,n_images);
        for kl=1:9
            for sign=[1, -1]
                % Gij = v(k)_pi*v(l)_pj.'
                Hij = R(:,:,est_rel_rots(p_i,p_j))*sign*e_kl(:,:,kl)*R(:,:,est_rel_rots(p_j,p_i)).';
                if u_G(pair_ind) < 0
                    Hij = J*Hij*J;
                end
                if norm(Hij_1-Hij) < min_norm_1
                    H(3*p_i-2:3*p_i,3*p_j-2:3*p_j,1) = Hij;
                    min_norm_1 = norm(Hij_1-Hij);
                end
                if norm(Hij_2-Hij) < min_norm_2
                    H(3*p_i-2:3*p_i,3*p_j-2:3*p_j,2) = Hij;
                    min_norm_2 = norm(Hij_2-Hij);
                end
                if norm(Hij_3-Hij) < min_norm_3
                    H(3*p_i-2:3*p_i,3*p_j-2:3*p_j,3) = Hij;
                    min_norm_3 = norm(Hij_3-Hij);
                end
            end
        end
    end
end

% blocks (i,i) i=1,..,n_images % 
for m=1:3
    for p_i=1:n_images
        for p_j=1:n_images
            % if p_i < p_j: H(p_i,p_j) = v(k)_pi*v(l)_pj.'
            %               H(p_i,p_i) = H(p_i,p_j)*H(p_i,p_j).' 
            % if p_i > p_j: H(p_j,p_i) = v(k)_pj*v(l)_pi.'
            %               H(p_i,p_i) = H(p_j,p_i).'*H(p_j,p_i) 
            if p_i < p_j 
                H(3*p_i-2:3*p_i,3*p_i-2:3*p_i,m) = ...
                    H(3*p_i-2:3*p_i,3*p_i-2:3*p_i,m) + ...
                    H(3*p_i-2:3*p_i,3*p_j-2:3*p_j,m)*H(3*p_i-2:3*p_i,3*p_j-2:3*p_j,m).';
            end
            if p_i > p_j 
                H(3*p_i-2:3*p_i,3*p_i-2:3*p_i,m) = ...
                    H(3*p_i-2:3*p_i,3*p_i-2:3*p_i,m) + ...
                    H(3*p_j-2:3*p_j,3*p_i-2:3*p_i,m).'*H(3*p_j-2:3*p_j,3*p_i-2:3*p_i,m);
            end
        end
        Hii = H(3*p_i-2:3*p_i,3*p_i-2:3*p_i,m)/(2*(n_images-1));
        % taking the best rank 1 approximation using SVD
        [u,s,v] = svd(Hii);
        [max_sv,argmax_sv] = max(diag(s));      % maximal singular value
        H(3*p_i-2:3*p_i,3*p_i-2:3*p_i,m) = max_sv*u(:,argmax_sv)*v(:,argmax_sv).';
    end
end


% blocks (j,i), i<j, equal (i,j).' %
for m=1:3
   H(:,:,m) = H(:,:,m) + H(:,:,m).';
end


%%% factorizing H(:,:,m) using SVD decomposition %%%

% In a clean setting each H(:,:,m), m=1,2,3, is of rank 1 whose eigenvector 
% is the concatenation of one of the rows of all rotation matrices
V = zeros(3*n_images,3);
for m=1:3
    [v, d] = eigs(H(:,:,m), 20, 'lm');
    [evals, ind] = sort(diag(d), 'descend');
    evect1 = v(:,ind(1));
    log_message('top 5 eigenvalues = %.2f %.2f %.2f %.2f %.2f\n',evals(1),evals(2),evals(3),evals(4),evals(5));
    for i=1:n_images
        vi = evect1(3*i-2:3*i);     % vi is one of the rows of Ri.
        vi = vi/norm(vi);           % each row should be of rank 1.
        V(3*i-2:3*i,m) = vi.';
    end
end

% reshaping V to the form (3,3,n_images):
rots = [reshape(V(:,1),3,1,n_images) reshape(V(:,2),3,1,n_images) reshape(V(:,3),3,1,n_images)];

% rots = O*gt_rots (gt_rots up to a group element), O is orthogonal matrix
% if det(O)=-1 then we multiply rots by -1.
if det(rots(:,:,1)) < 0
    rots = -rots;
end


