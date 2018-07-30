function vis = estimate_third_rows_c2(Rijs,Rijgs,conf,nImages,is_use_weights)
%
% Find the third row of each rotation matrix.
% 
% Input parameters:
%   vijs       A 3x3xnchoose2 array where each 3x3 slice holds the
%              third rows outer product of the corresponding pair of
%              matrices.
%   viis       A 3x3xn array where the i-th 3x3 slice holds the outer
%              product of the third row of Ri with itself
% Output parameters:
%   vis        A 3xnImages matrix whose i-th column equals the 
%              transpose of the third row of the rotation matrix Ri.

assert(nchoosek(nImages,2) == size(Rijs,3));
assert(size(Rijs,3) == size(Rijgs,3));

if is_use_weights
    W = conf + conf'; %+ eye(K);
    D = sum(W,2);
    W = diag(D)\W;  % normalize every row    
else
    W = ones(nImages,nImages);
end

% V is a 3nx3n matrix whose (i,j)-th block of size 3x3 holds 
% the outer product vij
V = zeros(3*nImages,3*nImages);
for i=1:nImages
    for j=i+1:nImages
        ind = uppertri_ijtoind(i,j,nImages);
        vij = (Rijs(:,:,ind)+Rijgs(:,:,ind))/2;
        V((i-1)*3+1:i*3,(j-1)*3+1:j*3) = vij*W(i,j);
    end
end

V = V + V.'; % since vij^{T} = vji

V = V + eye(3*nImages); % since on the diagonal there is Ri^{T}*Ri = I

[v, d] = eigs(V, 20, 'la');
[evals, ind] = sort(diag(d), 'descend');
% In a clean setting V is of rank 1 whose eigenvector is the concatenation
% of the third rows of all rotation matrices
evect1 = v(:,ind(1));

log_message('Third row estimation sync matrix: first 5 eigenvalues=[%.2f %.2f %.2f %.2f %.2f]',evals(1),evals(2),evals(3),evals(4),evals(5));

% figure,
% bar(evals);
% set(gca, 'FontSize', 20);
% title('evals of S_{33}');

vis = zeros(3,nImages);
for i=1:nImages
    vi = evect1((i-1)*3+1:i*3);
    vi = vi/norm(vi); % each row should be of rank-1
    vis(:,i) =  vi.';
end

end