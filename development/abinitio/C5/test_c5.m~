nImages = 50;
n_badImags = nImages*3;
n_theta = 360;
snr = 100000000;
max_shift  = 0;
shift_step = 1;
[projs,refq] = generate_c4_images(nImages,snr,65,'GAUSSIAN',max_shift,shift_step);
rots_gt = zeros(3,3,nImages);
for i=1:nImages
    rots_gt(:,:,i) = q_to_rot(refq(:,i))';
end

refq_bad = qrand(n_badImags);
rots_bad = zeros(3,3,n_badImags);
for i=1:n_badImags
    rots_bad(:,:,i) = q_to_rot(refq_bad(:,i))';
end



rots_all = cat(3,rots_gt,rots_bad);
nImages_all = size(rots_all,3);

for i=1:nImages_all
    for j=i+1:nImages_all

end




S_3n = zeros(3*nImages,3*nImages);

for i=1:nImages
    for j=i+1:nImages
        Ri = rots_gt(:,:,i);
        Rj = rots_gt(:,:,j);
        S_3n(3*(i-1)+1:3*i,3*(j-1)+1:3*j) = Ri.'*Rj;
    end
end

S_3n = S_3n + S_3n.'+ eye(3*nImages);

[v, d] = eigs(S_3n, 20, 'la');
[evals, ind] = sort(diag(d), 'descend');
evects = v(:,ind(1:3));

disp('-first 10 eigenvalues---');
disp(evals(1:10));
disp('---------')

figure;
bar(evals);
set(gca,'FontSize', 20);
title('S3n');

TOL = 1e-12;
rotations = zeros(3,3,nImages);
for i = 1:nImages
    Rii = evects((i-1)*3+1:i*3,:);
    n2 = sum(Rii.^2, 2);
    Rii = Rii./repmat(sqrt(n2),1,3);
    
    [U,~,V]=svd(Rii); %approx by an O(3) matrix
    Rii=U*V.';
    
    if norm(Rii.'*Rii -eye(3), 'fro') > TOL || abs(abs(det(Rii)) -1) > TOL
        error('aaa');
    end
    rotations(:,:,i) =  Rii.';
end

[rot_alligned,err_in_degrees,mse] = analyze_results(rotations,n_theta,refq_gt);