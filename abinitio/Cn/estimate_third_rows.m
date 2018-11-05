function vis = estimate_third_rows(vijs,viis,n_symm)
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

if ~exist('is_conjugate_with_vii','var')
    % an empirical observation. But better check the alternative if reconstruction is not good enough
    if(n_symm==3 || n_symm==4)
        log_message('Not conjugating each vij with vii and vjj. Consider testing alternative if reconstruction is not satisfactory');
        is_conjugate_with_vii = false;
    else
        log_message('Conjugating each vij with vii and vjj. Consider testing alternative if reconstruction is not satisfactory');
        is_conjugate_with_vii = true;
    end
end

[nr,nc,nImages] = size(viis);
assert(nr == 3 && nc == 3);


[nr_,nc_,n_] = size(vijs);
assert(nr_ == 3 && nc_ == 3);

assert(n_ == nchoosek(nImages,2));

% V is a 3nx3n matrix whose (i,j)-th block of size 3x3 holds 
% the outer product vij
V = zeros(3*nImages,3*nImages);
for i=1:nImages
%     V((i-1)*3+1:i*3,(i-1)*3+1:i*3) = viis(:,:,i);
    for j=i+1:nImages
        ind = uppertri_ijtoind(i,j,nImages);
        if is_conjugate_with_vii
            V((i-1)*3+1:i*3,(j-1)*3+1:j*3) = viis(:,:,i)*vijs(:,:,ind)*viis(:,:,j);
        else
            V((i-1)*3+1:i*3,(j-1)*3+1:j*3) = vijs(:,:,ind);
        end
        
    end
end

V = V + V.'; % since vij^{T} = vji

for i=1:nImages
    V((i-1)*3+1:i*3,(i-1)*3+1:i*3) = viis(:,:,i);
end

[v, d] = eigs(V, 20, 'lm');
[evals, ind] = sort(diag(d), 'descend');
% In a clean setting V is of rank 1 whose eigenvector is the concatenation
% of the third rows of all rotation matrices
evect1 = v(:,ind(1));

log_message('first 5 eigenvalues=[%.2f %.2f %.2f %.2f %.2f]',evals(1),evals(2),evals(3),evals(4),evals(5));

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