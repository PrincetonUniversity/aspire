n = 13;
m = 32;

fprintf('------\n| 1D |\n------\n');

c = rand(n, 1) + i*rand(n, 1);
om = m*(rand(n, 1)-0.5);

tt = tic;
f1 = nufft_1d(c, om, 'double', m);
t1 = toc(tt);

tt = tic;
f2 = n*nufft1d1(n, 2*pi/m*om, c, 1, 1e-10, m);
t2 = toc(tt);

tt = tic;
f3 = slow_nufft_1d(c, om, m);
t3 = toc(tt);

fprintf('nufft_1d: %g s\n', t1);
fprintf('nufft1d1: %g s\n', t2);
fprintf('slow_nufft_1d: %g s\n', t3);

fprintf('error of nufft_1d: %g\n', norm(f1(:)-f3(:))/norm(f3(:)));
fprintf('error of nufft1d1: %g\n', norm(f2(:)-f3(:))/norm(f3(:)));

cf = rand(m, 1) + i*rand(m, 1);

tt = tic;
g1 = nufft_t_1d(cf, 2*pi/m*om, 'double');
t1 = toc(tt);

tt = tic;
g2 = nufft1d2(n, 2*pi/m*om, 1, 1e-10, m, cf);
t2 = toc(tt);

tt = tic;
g3 = slow_nufft_t_1d(cf, 2*pi/m*om);
t3 = toc(tt);

fprintf('nufft_1d_t: %g s\n', t1);
fprintf('nufft1d2: %g s\n', t2);
fprintf('slow_nufft_1d_t: %g s\n', t3);

fprintf('error of nufft_t_1d: %g\n', norm(g1(:)-g3(:))/norm(g3(:)));
fprintf('error of nufft1d2: %g\n', norm(g2(:)-g3(:))/norm(g3(:)));

fprintf('\n');

fprintf('------\n| 2D |\n------\n');

om = m*(rand(n, 2)-0.5);

tt = tic;
f1 = nufft_2d(c, om, 'double', m);
t1 = toc(tt);

tt = tic;
f2 = n*nufft2d1(n, 2*pi/m*om(:,1), 2*pi/m*om(:,2), c, 1, 1e-10, m, m);
t2 = toc(tt);

tt = tic;
f3 = slow_nufft_2d(c, om, m);
t3 = toc(tt);

fprintf('nufft_2d: %g s\n', t1);
fprintf('nufft2d1: %g s\n', t2);
fprintf('slow_nufft_2d: %g s\n', t3);

fprintf('error of nufft_2d: %g\n', norm(f1(:)-f3(:))/norm(f3(:)));
fprintf('error of nufft2d1: %g\n', norm(f2(:)-f3(:))/norm(f3(:)));

cf = rand(m*ones(1, 2)) + i*rand(m*ones(1, 2));

tt = tic;
g1 = nufft_t_2d(cf, 2*pi/m*om, 'double');
t1 = toc(tt);

tt = tic;
g2 = nufft2d2(n, 2*pi/m*om(:,1), 2*pi/m*om(:,2), 1, 1e-10, m, m, cf);
t2 = toc(tt);

tt = tic;
g3 = slow_nufft_t_2d(cf, 2*pi/m*om);
t3 = toc(tt);

tt = tic;
precomp = nufft_t_2d_prepare(2*pi/m*om, m, 'double');
t4p = toc(tt);
tt = tic;
g4 = nufft_t_2d_execute(cf, precomp);
t4 = toc(tt);

tt = tic;
g5 = nufft_t_2d_optimized(cf, 2*pi/m*om, 'double');
t5 = toc(tt);

tt = tic;
g6 = nufft_t_2d_noQ(cf, 2*pi/m*om, 'double');
t6 = toc(tt);

fprintf('nufft_2d_t: %g s\n', t1);
fprintf('nufft2d2: %g s\n', t2);
fprintf('slow_nufft_2d_t: %g s\n', t3);
fprintf('nufft_t_2d_execute: %g s (%g s)\n', t4p+t4, t4);
fprintf('nufft_t_2d_optimized: %g s\n', t5);
fprintf('nufft_t_2d_noQ: %g s\n', t6);

fprintf('error of nufft_t_2d: %g\n', norm(g1(:)-g3(:))/norm(g3(:)));
fprintf('error of nufft2d2: %g\n', norm(g2(:)-g3(:))/norm(g3(:)));
fprintf('error of nufft_t_2d_execute: %g\n', norm(g4(:)-g3(:))/norm(g3(:)));
fprintf('error of nufft_t_2d_optimized: %g\n', norm(g5(:)-g3(:))/norm(g3(:)));
fprintf('error of nufft_t_2d_noQ: %g\n', norm(g6(:)-g3(:))/norm(g3(:)));

fprintf('\n');

fprintf('------\n| 3D |\n------\n');

om = m*(rand(n, 3)-0.5);

tt = tic;
f1 = nufft_3d_ref(c, om, 'double', m);
t1 = toc(tt);

tt = tic;
f2 = n*nufft3d1(n, 2*pi/m*om(:,1), 2*pi/m*om(:,2), 2*pi/m*om(:,3), c, 1, 1e-10, m, m, m);
f2 = reshape(f2, m*ones(1, 3));
t2 = toc(tt);

tt = tic;
f3 = slow_nufft_3d(c, om, m);
t3 = toc(tt);

tt = tic;
f4 = nufft_3d(c, om, 'double', m);
t4 = toc(tt);

fprintf('nufft_3d_ref: %g s\n', t1);
fprintf('nufft3d1: %g s\n', t2);
fprintf('slow_nufft_3d: %g s\n', t3);
fprintf('nufft_3d: %g s\n', t4);

fprintf('error of nufft_3d_ref: %g\n', norm(f1(:)-f3(:))/norm(f3(:)));
fprintf('error of nufft3d1: %g\n', norm(f2(:)-f3(:))/norm(f3(:)));
fprintf('error of nufft_3d: %g\n', norm(f4(:)-f3(:))/norm(f3(:)));

cf = rand(m*ones(1, 3)) + i*rand(m*ones(1, 3));

tt = tic;
g1 = nufft_t_3d(cf, 2*pi/m*om, 'double');
t1 = toc(tt);

tt = tic;
g2 = nufft3d2(n, 2*pi/m*om(:,1), 2*pi/m*om(:,2), 2*pi/m*om(:,3), 1, 1e-10, m, m, m, cf);
t2 = toc(tt);

tt = tic;
g3 = slow_nufft_t_3d(cf, 2*pi/m*om);
t3 = toc(tt);

tt = tic;
precomp = nufft_t_3d_prepare_1(2*pi/m*om, m, 'double');
t4p = toc(tt);
tt = tic;
g4 = nufft_t_3d_execute_1(cf, precomp);
t4 = toc(tt);

fprintf('nufft_3d_t: %g s\n', t1);
fprintf('nufft3d2: %g s\n', t2);
fprintf('slow_nufft_3d_t: %g s\n', t3);
fprintf('nufft_t_3d_execute_1: %g s (%g s)\n', t4p+t4, t4);

fprintf('error of nufft_t_3d: %g\n', norm(g1(:)-g3(:))/norm(g3(:)));
fprintf('error of nufft3d2: %g\n', norm(g2(:)-g3(:))/norm(g3(:)));
fprintf('error of nufft_t_3d_execute_1: %g\n', norm(g4(:)-g3(:))/norm(g3(:)));

