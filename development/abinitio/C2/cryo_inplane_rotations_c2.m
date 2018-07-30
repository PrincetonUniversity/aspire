function rots = cryo_inplane_rotations_c2(vis,Rijs,Rijgs,is_use_weight,conf)

nImages = size(vis,2);

assert(size(Rijs,3) == nchoosek(nImages,2));
assert(size(Rijgs,3) == nchoosek(nImages,2));

if is_use_weight
    W = conf + conf' ; %+ eye(K);
    D = sum(W, 2);
    W = diag(D)\W;  % normalize every row
    % adding 'I' for debug purposes. Namely that in clean setting thw spectrum is indeed (1,0,0,...0). 
    %Otherwise everything is shifted by 1
    W = W + eye(nImages); 
else
    W = ones(nImages,nImages);
end


H = zeros(nImages,nImages);

Ris_tilde = zeros(3,3,nImages);
for i=1:nImages
    v_i = vis(:,i).';
    Ris_tilde(:,:,i) = complete_3rdRow_to_rot(v_i);
end

for i = 1:nImages
    
    Ri_tilde = Ris_tilde(:,:,i);
    
    for j = i+1:nImages
        
        Rj_tilde  = Ris_tilde(:,:,j);
        
        ind = uppertri_ijtoind(i,j,nImages);
        
        Rij  = Rijs(:,:,ind);
        Rijg = Rijgs(:,:,ind);
        
        Uij = Ri_tilde*Rij*Rj_tilde.'; %U_ij is x-y-plane rotation
        [u, ~, v] = svd(Uij(1:2, 1:2));
        Uij = u*v.';
        
        Uijg = Ri_tilde*Rijg*Rj_tilde.';
        [u, ~, v] = svd(Uijg(1:2, 1:2));
        Uijg = u*v.';
        
        U = (Uij*Uij + Uijg*Uijg)/2;
        [u, ~, v] = svd(U);
        U = u*v.';
        
        H(i,j) = U(1,1) + sqrt(-1)*U(2,1);
        H(i,j) = H(i,j)*W(i,j);
    end
end

% note entry (i,j) corresponds to exp^(-i(-theta_i+theta_i)). Therefore to
% construct the hermitian matrix, entry (j,i) is the **conjugate** of entry
% (i,j)
% transpose
H = H + H' + eye(nImages); % put 1 on diagonal since : exp^(-i*0) = 1

[v, d] = eigs(H, 10, 'lm');
[evals, ind] = sort(diag(d), 'descend');
evect1 = v(:,ind(1)); %zi = exp(-i*2*ai), Ri = Rot(ai)*Ri0

log_message('In-plane rotation sync matrix : first 5 eigenvalues=[%.2f %.2f %.2f %.2f %.2f]',evals(1),evals(2),evals(3),evals(4),evals(5));

% figure,
% hist(evals, 40);
% set(gca, 'FontSize', 20);
% title('SO(2) Sync');

rots = zeros(3,3,nImages);
for i  = 1:nImages
    zi  = evect1(i);
    zi  = zi/abs(zi); % rescale so it lies on unit circle
    c   = real(sqrt(zi));  % we have a double angle, so apply sqrt
    s   = -imag(sqrt(zi)); % we have a double angle, so apply sqrt.
    Rai = [c -s 0; s c 0; 0 0 1];
    
   Ri_tilde = Ris_tilde(:,:,i);
    
    rots(:,:,i) = Rai*Ri_tilde;
end

end

function R = complete_3rdRow_to_rot(r3)

% build R \in SO(3) from the 3rd Row
% 3rd col is not aligning to z-axis

tmp = sqrt(r3(1)^2 + r3(2)^2);

r1 = [r3(2)/tmp -r3(1)/tmp 0];
r2 = [r3(1)*r3(3)/tmp r3(2)*r3(3)/tmp -tmp];

R = [r1; r2; r3];

end