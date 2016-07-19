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

fprintf('error of nufft_1d vs nufft1d1: %g\n', norm(f1(:)-f2(:))/norm(f2(:)));
fprintf('error of nufft_1d vs slow_nufft_1d: %g\n', norm(f1(:)-f3(:))/norm(f3(:)));

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

fprintf('error of nufft_t_1d vs nufft1d2: %g\n', norm(g1(:)-g2(:))/norm(g2(:)));
fprintf('error of nufft_t_1d vs slow_nufft_t_1d: %g\n', norm(g1(:)-g3(:))/norm(g3(:)));

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

fprintf('error of nufft_2d vs nufft2d1: %g\n', norm(f1(:)-f2(:))/norm(f2(:)));
fprintf('error of nufft_2d vs slow_nufft_2d: %g\n', norm(f1(:)-f3(:))/norm(f3(:)));

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

fprintf('nufft_2d_t: %g s\n', t1);
fprintf('nufft2d2: %g s\n', t2);
fprintf('slow_nufft_2d_t: %g s\n', t3);

fprintf('error of nufft_t_2d vs nufft2d2: %g\n', norm(g1(:)-g2(:))/norm(g2(:)));
fprintf('error of nufft_t_2d vs slow_nufft_t_2d: %g\n', norm(g1(:)-g3(:))/norm(g3(:)));

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

fprintf('nufft_3d: %g s\n', t1);
fprintf('nufft3d1: %g s\n', t2);
fprintf('slow_nufft_3d: %g s\n', t3);

fprintf('error of nufft_3d vs nufft3d1: %g\n', norm(f1(:)-f2(:))/norm(f2(:)));
fprintf('error of nufft_3d vs slow_nufft_3d: %g\n', norm(f1(:)-f3(:))/norm(f3(:)));

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

fprintf('nufft_3d_t: %g s\n', t1);
fprintf('nufft3d2: %g s\n', t2);
fprintf('slow_nufft_3d_t: %g s\n', t3);

fprintf('error of nufft_t_3d vs nufft3d2: %g\n', norm(g1(:)-g2(:))/norm(g2(:)));
fprintf('error of nufft_t_3d vs slow_nufft_t_3d: %g\n', norm(g1(:)-g3(:))/norm(g3(:)));

