projs=ReadMRC('~/tmp/80S_89/averages_nn50_group1.mrc');
projs=projs(:,:,1:500);

open_log(0);

n_r=floor(size(projs,1)/2);
n_theta=360;
[pf,~]=cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections     

N=size(projs,3);
pf=cryo_raynormalize(reshape(pf,n_r,n_theta*N));
pf=reshape(pf,n_r,n_theta,N);


psd=cryo_noise_estimation_pfs(projs,n_r,n_theta);

M=ceil(1/(2*(sind(1/4))^2));
max_itr=4;
verbose=1;
KNN=50;
max_shift_1d=10;
shift_step_1d=1;
tic;
[clstack_ref]=cryo_clmatrix_ml_ref(pf,M,KNN,diag(psd),max_itr,verbose,max_shift_1d,shift_step_1d);
t_ref=toc
tic;
[clstack]=cryo_clmatrix_ml(pf,M,KNN,diag(psd),max_itr,verbose,max_shift_1d,shift_step_1d);
t_ml=toc
tic;
[clstack_gpu]=cryo_clmatrix_ml_gpu(pf,M,KNN,diag(psd),max_itr,verbose,max_shift_1d,shift_step_1d);
t_gpu=toc
