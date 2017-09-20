% test the function cryo_clmatrix_gpu
%
% Compare cryo_clmatrix to cryo_clmatrix_gpu bit by bit, when applying
% averaging to the correlation map.
%
% Yoel Shkolnisky, June 2015.
%

K=25;  % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=3;     % The shifts in each projection are between -max_shift_2d  
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
shift_step_1d=1;    % 1D shifts steps to use in common line matching.
thetaTol=5;  % A common line is considered correct if it deviates from the 
             % true common line by up to thetaTol degrees.
snr=1;  % Not used. Only clean proejctions are used.
filter_radius=5; % How many nearby correlation values to average.

% Generate proejctions. Use only the generated clean projections "projs".
initstate;
[projs,~,~,rots_ref]=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

% take Fourier transform of projections.
n_theta=360;  
n_r=ceil(n/2);
pf=cryo_pft(projs,n_r,n_theta,'single');  

max_shift_1d=ceil(2*sqrt(2)*max_shift_2d); % Maximal 1D shift that may be 
    % introduced by a 2D shift of max_shift_2d.
open_log(0);

% Find common lines.
tic;
% Switch GPU to DOUBLE for percise match of CPU and GPU versions.
% But then the GPU code is slower.
[clstack_gpu,corrstack_gpu,shift_equations_gpu,shift_equations_map_gpu]...
    = cryo_clmatrix_gpu(pf,K,1,max_shift_1d,shift_step_1d,filter_radius);
t_gpu=toc;
tic;
[clstack_cpu,corrstack_cpu,shift_equations_cpu,shift_equations_map_cpu] ...
    = cryo_clmatrix(pf,K,1,max_shift_1d,shift_step_1d,filter_radius);
t_cpu=toc;

% Check all output variables.
% All errors should be of the order of machine precision used for GPU
% calculations.
fprintf('clstack error = %7.4e\n',relerr(clstack_gpu,clstack_cpu));
fprintf('corrstack error = %7.4e\n',relerr(corrstack_gpu,corrstack_cpu));
shift_diff=shift_equations_gpu-shift_equations_cpu;
err=full(sqrt(sum(shift_diff(:).^2)));
fprintf('shift_equations error = %7.4e\n',err);
shift_diff=shift_equations_map_gpu-shift_equations_map_cpu;
err=full(sqrt(sum(shift_diff(:).^2)));
fprintf('shift_equations_map error = %7.4e\n',err);
fprintf('Timing t_cpu/t_gpu = %7.4e\n',t_cpu/t_gpu);
