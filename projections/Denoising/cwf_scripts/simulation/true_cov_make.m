
cd ~/aspire
initpath

cd ~/Dropbox/asinger/80S/
initpath_t
addpath ~/aspire/io/
addpath ~/aspire/io/log/
addpath ~/aspire/common/
addpath ~/aspire/projections/simulation
addpath ~/aspire/projections/
addpath ~/aspire/projections/class_average/preprocess/
addpath ~/aspire/projections/epsd/  

c=0.48; R=45; % Same as in fred_shrinkage experiment, for consistency in ground truth
n_r = ceil(4*c*R);
tic_basis=tic;
[ basis, sample_points ] = precomp_fb( n_r, R, c );
timing.basis=toc(tic_basis)
num_pool=20;

% Calculate mean
tr_mean=zeros(105,105);
for i=1:10
	tic_mean=tic;
	filename=fullfile('/scratch/tbhamre/1e6_proj_6454/', sprintf('set%d',i));
	load (filename);
	disp('Loaded')
	proj=mask_corners(proj);
	hatI_curr=cfft2(proj);
	tr_mean=tr_mean+mean(hatI_curr,3);	
	clear proj
	toc_mean=toc(tic_mean)
end
tr_mean=tr_mean/10;
% calculate covariance

tr_cov=cell(length(unique(basis.ang_freqs)), 1);

for i=1:10
	tic_cov=tic;
	filename=fullfile('/scratch/tbhamre/1e6_proj_6454/', sprintf('set%d',i));
	load (filename);
	sprintf('Loaded %d',i)
	N_images=100000;
	proj=mask_corners(proj);
	hatI_curr=cfft2(proj);
	nim=1e6;
	L0=size(hatI_curr,1);
	[ tr_coeff_pos_k ] = coeff_demean( icfft2(bsxfun(@minus,hatI_curr,tr_mean)) , R, basis, sample_points, num_pool);
	for k=unique(basis.ang_freqs)'
		if k>=0
			tr_c=tr_coeff_pos_k{k+1};
			if (i==1)
				tr_cov{k+1}=real(tr_c*tr_c')/(nim);
			else
				tr_cov{k+1}=tr_cov{k+1}+real(tr_c*tr_c')/(nim);
			end
	   	end
	end
	toc_cov=toc(tic_cov)
end
%% Save results
fpath = sprintf('/scratch/tbhamre/jsb_expt/sim/tr_cov_1e6.mat');
save(fpath,'tr_cov', 'tr_coeff_pos_k','-v7.3');


