% Test the function cryo_add shifts.
%
% Yoel Shkolnisky, September 2013.

K=10;
SNR=1;
max_shift=5;
shift_step=1;
[clean_projections_1,noisy_projections_1,ref_shifts,ref_q]=gen_projections_v2(K,SNR,max_shift,shift_step);
n=size(clean_projections_1,1);

tic;
[clean_projections_2, noisy_projections_2, shifts, q] = ...
    cryo_gen_projections(n,K,SNR,[],[],ref_shifts,ref_q,'single');
toc
figure(1);
subplot(2,3,1); imagesc(clean_projections_1(:,:,1)); colorbar; axis image;
subplot(2,3,2); imagesc(clean_projections_2(:,:,1)); colorbar; axis image;
subplot(2,3,3); imagesc(clean_projections_1(:,:,1)-clean_projections_2(:,:,1)); colorbar; axis image;
err=norm(clean_projections_1(:)-clean_projections_2(:))/norm(clean_projections_1(:));
disp(err); % Should be of the order of machine precision 
           % (single or double machine precision according to the precision
           % passed to cryo_gen_projections above). 

subplot(2,3,4); imagesc(noisy_projections_1(:,:,1)); colorbar; axis image;
subplot(2,3,5); imagesc(noisy_projections_2(:,:,1)); colorbar; axis image;
subplot(2,3,6); imagesc(noisy_projections_1(:,:,1)-noisy_projections_2(:,:,1)); colorbar; axis image;
err=norm(noisy_projections_1(:)-noisy_projections_2(:))/norm(noisy_projections_1(:));
disp(err); % The result is not be small since the two functions define the 
           % SNR differently. To get same results up to machine precision,
           % uncomment lines 54-55 in cryo_addnoise.
           % The figure will always show tiny error since it is used to
           % determine the SNR. See cryo_addnoise for details.
