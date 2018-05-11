% Test the function cryo_masking_radius_2d.
% Generate K images, and plot a bounding circle of the molecule in each
% image.
% Press any key to continue after each image.
%
% Yoel Shkolnisky, October 2016.

n=89;
K=10;
SNR=1;
initstate;
[~,noisy_projs]=cryo_gen_projections(n,K,SNR);

for k=1:K
    fprintf('Testing image %2d...',k);
    P=noisy_projs(:,:,k);
    rmin=cryo_masking_radius_2d(P,1);
    fprintf('r=%d\n',rmin);
    pause
end

