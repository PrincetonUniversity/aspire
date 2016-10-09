% Test the function cryo_masking_radius_3d.
%
% Yoel Shkolnisky, October 2016.

%mapname=cryo_fetch_emdID(2275);
%mapname=cryo_fetch_emdID(6487);
%map=ReadMRC(mapname);
map=single(map);
map=cryo_downsample(map,[89,89,89]);
SNR=10;
sigma=sqrt(var(map(:))/SNR);
noise=randn(size(map));
noisy_map=map+sigma*noise;
rmin=cryo_masking_radius_3d(noisy_map,1);
