N = 10;
M = 87;

old_libs = get_nufft_libraries();

nudft_warning = warning('query', 'aspire:using_nudft');
warning('off', nudft_warning.identifier);

fprintf('Testing NUFFT wrappers\n');

fprintf('N = %d, M = %d\n', N, M);

vol = randn(N*ones(1, 3))+i*randn(N*ones(1, 3));
fourier_pts = 2*pi*(rand(3, M)-0.5);

set_nufft_libraries('chemnitz');
tt = tic;
vol1_f = nufft3(vol, fourier_pts);
tm1_3d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft3 ''chemnitz'':', tm1_3d);

set_nufft_libraries('cims');
tt = tic;
vol2_f = nufft3(vol, fourier_pts);
tm2_3d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft3 ''cims'':', tm2_3d);

set_nufft_libraries('dft');
tt = tic;
vol3_f = nufft3(vol, fourier_pts);
tm3_3d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft3 ''dft'':', tm3_3d);

set_nufft_libraries('finufft');
tt = tic;
vol4_f = nufft3(vol, fourier_pts);
tm4_3d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft3 ''finufft'':', tm4_3d);

tt = tic;
vol0_f = nudft3(vol, fourier_pts);
tm0_3d = toc(tt);
fprintf('%-40s%15f s\n', 'nudft3:', tm0_3d);

err1_3d = norm(vol1_f(:)-vol0_f(:))/norm(vol0_f(:));
err2_3d = norm(vol2_f(:)-vol0_f(:))/norm(vol0_f(:));
err3_3d = norm(vol3_f(:)-vol0_f(:))/norm(vol0_f(:));
err4_3d = norm(vol4_f(:)-vol0_f(:))/norm(vol0_f(:));

fprintf('%-40s%15g\n', 'error nufft3 ''chemnitz'':', err1_3d);
fprintf('%-40s%15g\n', 'error nufft3 ''cims'':', err2_3d);
fprintf('%-40s%15g\n', 'error nufft3 ''dft'':', err3_3d);
fprintf('%-40s%15g\n', 'error nufft3 ''finufft'':', err4_3d);

vol_f = randn(M, 1)+i*randn(M, 1);

set_nufft_libraries('chemnitz');
tt = tic;
vol1 = anufft3(vol_f, fourier_pts, N*ones(1, 3));
tm1_a3d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft3 ''chemnitz'':', tm1_a3d);

set_nufft_libraries('cims');
tt = tic;
vol2 = anufft3(vol_f, fourier_pts, N*ones(1, 3));
tm2_a3d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft3 ''cims'':', tm2_a3d);

set_nufft_libraries('dft');
tt = tic;
vol3 = anufft3(vol_f, fourier_pts, N*ones(1, 3));
tm3_a3d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft3 ''dft'':', tm3_a3d);

set_nufft_libraries('finufft');
tt = tic;
vol4 = anufft3(vol_f, fourier_pts, N*ones(1, 3));
tm4_a3d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft3 ''finufft'':', tm4_a3d);

tt = tic;
vol0 = anudft3(vol_f, fourier_pts, N*ones(1, 3));
tm0_a3d = toc(tt);
fprintf('%-40s%15f s\n', 'anudft3:', tm0_a3d);

err1_a3d = norm(vol1(:)-vol0(:))/norm(vol0(:));
err2_a3d = norm(vol2(:)-vol0(:))/norm(vol0(:));
err3_a3d = norm(vol3(:)-vol0(:))/norm(vol0(:));
err4_a3d = norm(vol4(:)-vol0(:))/norm(vol0(:));

fprintf('%-40s%15g\n', 'error anufft3 ''chemnitz'':', err1_a3d);
fprintf('%-40s%15g\n', 'error anufft3 ''cims'':', err2_a3d);
fprintf('%-40s%15g\n', 'error anufft3 ''dft'':', err3_a3d);
fprintf('%-40s%15g\n', 'error anufft3 ''finufft'':', err4_a3d);

im = randn(N*ones(1, 2))+i*randn(N*ones(1, 2));
fourier_pts = 2*pi*(rand(2, M)-0.5);

set_nufft_libraries('chemnitz');
tt = tic;
im1_f = nufft2(im, fourier_pts);
tm1_2d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft2 ''chemnitz'':', tm1_2d);

set_nufft_libraries('cims');
tt = tic;
im2_f = nufft2(im, fourier_pts);
tm2_2d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft2 ''cims'':', tm2_2d);

set_nufft_libraries('dft');
tt = tic;
im3_f = nufft2(im, fourier_pts);
tm3_2d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft2 ''dft'':', tm3_2d);

set_nufft_libraries('finufft');
tt = tic;
im4_f = nufft2(im, fourier_pts);
tm4_2d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft2 ''finufft'':', tm4_2d);

tt = tic;
im0_f = nudft2(im, fourier_pts);
tm0_2d = toc(tt);
fprintf('%-40s%15f s\n', 'nudft2:', tm0_2d);

err1_2d = norm(im1_f(:)-im0_f(:))/norm(im0_f(:));
err2_2d = norm(im2_f(:)-im0_f(:))/norm(im0_f(:));
err3_2d = norm(im3_f(:)-im0_f(:))/norm(im0_f(:));
err4_2d = norm(im4_f(:)-im0_f(:))/norm(im0_f(:));

fprintf('%-40s%15g\n', 'error nufft2 ''chemnitz'':', err1_2d);
fprintf('%-40s%15g\n', 'error nufft2 ''cims'':', err2_2d);
fprintf('%-40s%15g\n', 'error nufft2 ''dft'':', err3_2d);
fprintf('%-40s%15g\n', 'error nufft2 ''finufft'':', err4_2d);

im_f = randn(M, 1)+i*randn(M, 1);

set_nufft_libraries('chemnitz');
tt = tic;
im1 = anufft2(im_f, fourier_pts, N*ones(1, 2));
tm1_a2d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft2 ''chemnitz'':', tm1_a2d);

set_nufft_libraries('cims');
tt = tic;
im2 = anufft2(im_f, fourier_pts, N*ones(1, 2));
tm2_a2d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft2 ''cims'':', tm2_a2d);

set_nufft_libraries('dft');
tt = tic;
im3 = anufft2(im_f, fourier_pts, N*ones(1, 2));
tm3_a2d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft2 ''dft'':', tm3_a2d);

set_nufft_libraries('finufft');
tt = tic;
im4 = anufft2(im_f, fourier_pts, N*ones(1, 2));
tm4_a2d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft2 ''finufft'':', tm4_a2d);

tt = tic;
im0 = anudft2(im_f, fourier_pts, N*ones(1, 2));
tm0_a2d = toc(tt);
fprintf('%-40s%15f s\n', 'anudft2:', tm0_a2d);

err1_a2d = norm(im1(:)-im0(:))/norm(im0(:));
err2_a2d = norm(im2(:)-im0(:))/norm(im0(:));
err3_a2d = norm(im3(:)-im0(:))/norm(im0(:));
err4_a2d = norm(im4(:)-im0(:))/norm(im0(:));

fprintf('%-40s%15g\n', 'error anufft2 ''chemnitz'':', err1_a2d);
fprintf('%-40s%15g\n', 'error anufft2 ''cims'':', err2_a2d);
fprintf('%-40s%15g\n', 'error anufft2 ''dft'':', err3_a2d);
fprintf('%-40s%15g\n', 'error anufft2 ''finufft'':', err4_a2d);

sig = randn(N, 1)+i*randn(N, 1);
fourier_pts = 2*pi*(rand(1, M)-0.5);

set_nufft_libraries('chemnitz');
tt = tic;
sig1_f = nufft1(sig, fourier_pts);
tm1_1d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft1 ''chemnitz'':', tm1_1d);

set_nufft_libraries('cims');
tt = tic;
sig2_f = nufft1(sig, fourier_pts);
tm2_1d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft1 ''cims'':', tm2_1d);

set_nufft_libraries('dft');
tt = tic;
sig3_f = nufft1(sig, fourier_pts);
tm3_1d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft1 ''dft'':', tm3_1d);

set_nufft_libraries('dft');
tt = tic;
sig4_f = nufft1(sig, fourier_pts);
tm4_1d = toc(tt);
fprintf('%-40s%15f s\n', 'nufft1 ''finufft'':', tm4_1d);

tt = tic;
sig0_f = nudft1(sig, fourier_pts);
tm0_1d = toc(tt);
fprintf('%-40s%15f s\n', 'nudft1:', tm0_1d);

err1_1d = norm(sig1_f(:)-sig0_f(:))/norm(sig0_f(:));
err2_1d = norm(sig2_f(:)-sig0_f(:))/norm(sig0_f(:));
err3_1d = norm(sig3_f(:)-sig0_f(:))/norm(sig0_f(:));
err4_1d = norm(sig4_f(:)-sig0_f(:))/norm(sig0_f(:));

fprintf('%-40s%15g\n', 'error nufft1 ''chemnitz'':', err1_1d);
fprintf('%-40s%15g\n', 'error nufft1 ''cims'':', err2_1d);
fprintf('%-40s%15g\n', 'error nufft1 ''dft'':', err3_1d);
fprintf('%-40s%15g\n', 'error nufft1 ''finufft'':', err4_1d);

sig_f = randn(M, 1)+i*randn(M, 1);

set_nufft_libraries('chemnitz');
tt = tic;
sig1 = anufft1(sig_f, fourier_pts, N);
tm1_a1d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft1 ''chemnitz'':', tm1_a1d);

set_nufft_libraries('cims');
tt = tic;
sig2 = anufft1(sig_f, fourier_pts, N);
tm2_a1d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft1 ''cims'':', tm2_a1d);

set_nufft_libraries('dft');
tt = tic;
sig3 = anufft1(sig_f, fourier_pts, N);
tm3_a1d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft1 ''dft'':', tm3_a1d);

set_nufft_libraries('finufft');
tt = tic;
sig4 = anufft1(sig_f, fourier_pts, N);
tm4_a1d = toc(tt);
fprintf('%-40s%15f s\n', 'anufft1 ''finufft'':', tm4_a1d);

tt = tic;
sig0 = anudft1(sig_f, fourier_pts, N);
tm0_a1d = toc(tt);
fprintf('%-40s%15f s\n', 'anudft1:', tm0_a1d);

err1_a1d = norm(sig1(:)-sig0(:))/norm(sig0(:));
err2_a1d = norm(sig2(:)-sig0(:))/norm(sig0(:));
err3_a1d = norm(sig3(:)-sig0(:))/norm(sig0(:));
err4_a1d = norm(sig4(:)-sig0(:))/norm(sig0(:));

fprintf('%-40s%15g\n', 'error anufft1 ''chemnitz'':', err1_a1d);
fprintf('%-40s%15g\n', 'error anufft1 ''cims'':', err2_a1d);
fprintf('%-40s%15g\n', 'error anufft1 ''dft'':', err3_a1d);
fprintf('%-40s%15g\n', 'error anufft1 ''finufft'':', err4_a1d);

set_nufft_libraries(old_libs);
warning(nudft_warning.state, nudft_warning.identifier);
