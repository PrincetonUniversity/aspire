function [est_inv_rots] = deterministic_rounding(Gram)
% Deterministic rounding procedure to recover the rotations from the Gram
% matrix.
%
% Input:
%   Gram: a 2K x 2K Gram matrix.
%
% Output:
%   est_inv_rots:  a 3 x 3 x K matrix, the ith estimated inverse rotation
%       matrix is est_inv_rots(:,:,i).
%
% Lanhui Wang, Aug 8, 2013

K = size(Gram, 1)/2; % number of projections;
[V,D] = eigs(Gram,4);
eig_vals = diag(D);
[sorted_eigs, sort_idx] = sort(real(eig_vals),'descend');
disp('Four largest eigenvalues: ');
sorted_eigs(1:4)
VV(:,(1:4)) = V(:,sort_idx(1:4));
V1 = VV(1:K,1:3);
V2 = VV(K+1:2*K,1:3);

V1 = V1';
V2 = V2';
A=ATA_solver(V1,V2,K);
R1 = A*V1;
R2 = A*V2;

est_inv_rots = zeros(3,3,K);

for k=1:K;
    R = [R1(:,k) R2(:,k)];
    [U,~,V] = svd(R,0);
    R = U*V';
    est_inv_rots(:,:,k) = [R(:,1) R(:,2) cross(R(:,1),R(:,2))];
end;