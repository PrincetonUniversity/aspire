function test_pft
%
% Test the function pft.
%
% Yoel Shkolnisky, December 2009.

n=64;
n_r=50;
n_theta=50;

im=phantom(n);

% Use slow_ptf
tic;
pf1=slow_pft(im,n_r,n_theta);
t1=toc;

% Use pft
tic;
pf2=pft(im,n_r,n_theta,'single');
t2=toc;

err=norm(pf1-pf2)/norm(pf1);
fprintf('err=%8.6e  t1=%7.5fs  t2=%7.5fs\n',err,t1,t2);

