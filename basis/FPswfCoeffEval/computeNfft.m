function [ images_nufft ] = computeNfft( images, usFftPts, L, points_inside_the_circle)
% This function computes the NFFT of the images at the designated nodes

images_full = zeros((2*L+1)^2, size(images, 2));
images_full(points_inside_the_circle,:) = images;
images_full = vec_to_im(images_full);

% NOTE: For backwards compatibility, we have to transpose the x and y axes
% here, but this is probably not the right thing to do.
images_nufft = nufft2(permute(images_full, [2 1 3]), usFftPts.');
end

