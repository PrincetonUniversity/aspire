N = 10;
M = 87;

fprintf('Testing NUFFT wrappers\n');

fprintf('N = %d, M = %d\n', N, M);

vol = randn(N*ones(1, 3))+i*randn(N*ones(1, 3));
fourier_pts = 2*pi*(rand(3, M)-0.5);

set_nufft_libraries('chemnitz');
tt = tic;
vol1_f = nufft3(vol, fourier_pts);
tm1_3d = toc(tt);
fprintf('nufft3 ''chemnitz'': %f s\n', tm1_3d);

set_nufft_libraries('cims');
tt = tic;
vol2_f = nufft3(vol, fourier_pts);
tm2_3d = toc(tt);
fprintf('nufft3 ''cims'': %f s\n', tm2_3d);

set_nufft_libraries('dft');
tt = tic;
vol3_f = nufft3(vol, fourier_pts);
tm3_3d = toc(tt);
fprintf('nufft3 ''dft'': %f s\n', tm3_3d);

err1_3d = norm(vol1_f(:)-vol3_f(:))/norm(vol3_f(:));
err2_3d = norm(vol2_f(:)-vol3_f(:))/norm(vol3_f(:));

fprintf('error nufft3 ''chemnitz'': %g\n', err1_3d);
fprintf('error nufft3 ''cims'': %g\n', err2_3d);

vol_f = randn(M, 1)+i*randn(M, 1);

set_nufft_libraries('chemnitz');
tt = tic;
vol1 = anufft3(vol_f, fourier_pts, N*ones(1, 3));
tm1_a3d = toc(tt);
fprintf('anufft3 ''chemnitz'': %f s\n', tm1_a3d);

set_nufft_libraries('cims');
tt = tic;
vol2 = anufft3(vol_f, fourier_pts, N*ones(1, 3));
tm2_a3d = toc(tt);
fprintf('anufft3 ''cims'': %f s\n', tm2_a3d);

set_nufft_libraries('dft');
tt = tic;
vol3 = anufft3(vol_f, fourier_pts, N*ones(1, 3));
tm3_a3d = toc(tt);
fprintf('anufft3 ''dft'': %f s\n', tm3_a3d);

err1_a3d = norm(vol1(:)-vol3(:))/norm(vol3(:));
err2_a3d = norm(vol2(:)-vol3(:))/norm(vol3(:));

fprintf('error anufft3 ''chemnitz'': %g\n', err1_a3d);
fprintf('error anufft3 ''cims'': %g\n', err2_a3d);

im = randn(N*ones(1, 2))+i*randn(N*ones(1, 2));
fourier_pts = 2*pi*(rand(2, M)-0.5);

set_nufft_libraries('chemnitz');
tt = tic;
im1_f = nufft2(im, fourier_pts);
tm1_2d = toc(tt);
fprintf('nufft2 ''chemnitz'': %f s\n', tm1_2d);

set_nufft_libraries('cims');
tt = tic;
im2_f = nufft2(im, fourier_pts);
tm2_2d = toc(tt);
fprintf('nufft2 ''cims'': %f s\n', tm2_2d);

set_nufft_libraries('dft');
tt = tic;
im3_f = nufft2(im, fourier_pts);
tm3_2d = toc(tt);
fprintf('nufft2 ''dft'': %f s\n', tm2_2d);

err1_2d = norm(im1_f(:)-im3_f(:))/norm(im3_f(:));
err2_2d = norm(im2_f(:)-im3_f(:))/norm(im3_f(:));

fprintf('error nufft2 ''chemnitz'': %g\n', err1_2d);
fprintf('error nufft2 ''cims'': %g\n', err2_2d);

im_f = randn(M, 1)+i*randn(M, 1);

set_nufft_libraries('chemnitz');
tt = tic;
im1 = anufft2(im_f, fourier_pts, N*ones(1, 2));
tm1_a2d = toc(tt);
fprintf('anufft2 ''chemnitz'': %f s\n', tm1_a2d);

set_nufft_libraries('cims');
tt = tic;
im2 = anufft2(im_f, fourier_pts, N*ones(1, 2));
tm1_a2d = toc(tt);
fprintf('anufft2 ''cims'': %f s\n', tm1_a2d);

set_nufft_libraries('dft');
tt = tic;
im3 = anufft2(im_f, fourier_pts, N*ones(1, 2));
tm1_a2d = toc(tt);
fprintf('anufft2 ''dft'': %f s\n', tm1_a2d);

err1_a2d = norm(im1(:)-im3(:))/norm(im3(:));
err2_a2d = norm(im2(:)-im3(:))/norm(im3(:));

fprintf('error anufft2 ''chemnitz'': %g\n', err1_a2d);
fprintf('error anufft2 ''cims'': %g\n', err2_a2d);

sig = randn(N, 1)+i*randn(N, 1);
fourier_pts = 2*pi*(rand(1, M)-0.5);

set_nufft_libraries('chemnitz');
tt = tic;
sig1_f = nufft1(sig, fourier_pts);
tm1_1d = toc(tt);
fprintf('nufft1 ''chemnitz'': %f s\n', tm1_1d);

set_nufft_libraries('cims');
tt = tic;
sig2_f = nufft1(sig, fourier_pts);
tm2_1d = toc(tt);
fprintf('nufft1 ''cims'': %f s\n', tm1_2d);

set_nufft_libraries('dft');
tt = tic;
sig3_f = nufft1(sig, fourier_pts);
tm3_1d = toc(tt);
fprintf('nufft1 ''dft'': %f s\n', tm1_3d);

err1_1d = norm(sig1_f(:)-sig3_f(:))/norm(sig3_f(:));
err2_1d = norm(sig2_f(:)-sig3_f(:))/norm(sig3_f(:));

fprintf('error nufft1 ''chemnitz'': %g\n', err1_1d);
fprintf('error nufft1 ''cims'': %g\n', err2_1d);

sig_f = randn(M, 1)+i*randn(M, 1);

set_nufft_libraries('chemnitz');
tt = tic;
sig1 = anufft1(sig_f, fourier_pts, N);
tm1_a1d = toc(tt);
fprintf('anufft1 ''chemnitz'': %f s\n', tm1_a1d);

set_nufft_libraries('cims');
tt = tic;
sig2 = anufft1(sig_f, fourier_pts, N);
tm2_a1d = toc(tt);
fprintf('anufft1 ''cims'': %f s\n', tm2_a1d);

set_nufft_libraries('dft');
tt = tic;
sig3 = anufft1(sig_f, fourier_pts, N);
tm3_a1d = toc(tt);
fprintf('anufft1 ''dft'': %f s\n', tm3_a1d);

err1_a1d = norm(sig1(:)-sig3(:))/norm(sig3(:));
err2_a1d = norm(sig2(:)-sig3(:))/norm(sig3(:));

fprintf('error anufft1 ''chemnitz'': %g\n', err1_a1d);
fprintf('error anufft1 ''cims'': %g\n', err2_a1d);
