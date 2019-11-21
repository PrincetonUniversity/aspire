function [ f ] = nufft_cleanrib( vol1, azi, polar, r )
% gives the fourier coefficients of clean ribosome 70S at sampling points
% of (azi, polar) column vectors
%el = pi/2 - polar;
%[az1,el1] = ndgrid(azi,el);
n = size(azi,1);
[x,y,z] = sph2cart(azi,polar,ones(n,n)*r);
carts = [x(:) y(:) z(:)];
rate = size(vol1,1);
Nd = [rate rate rate];
Jd = Nd/2;
Kd = Nd*2;
n_shift = Nd/2;
st = nufft_init(carts,Nd,Jd,Kd,n_shift,'table', 2^11, 'minmax:kb');
f = nufft(vol1,st);

end

