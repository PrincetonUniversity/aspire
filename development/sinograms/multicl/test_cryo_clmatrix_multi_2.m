% Compare the function cryo_clmatrix_multi to cryo_clmatrix under noise
%
% Generate noisy projections, find common lines using cryo_clmatrix_multi
% and cryo_clmatrix, and compare the results.
%
% Images are of odd size.
%
% Yoel Shkolnisky, October 2013.


%% Odd image, no shifts
K=10;  % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=0;     % The shifts in each projection are between -max_shift_2d  
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
shift_step_1d=1;    % 1D shifts steps to use in common line matching.
thetaTol=5;  % A common line is considered correct if it deviates from the 
             % true common line by up to thetaTol degrees.
snr=1/10;  % Not used. Only clean proejctions are used.

% Generate noisy proejctions. 
[~,noisy_projs,~,refq]=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

% take Fourier transform of projections.
n_theta=360;  
n_r=ceil(n/2);
pf=cryo_pft(noisy_projs,n_r,n_theta,'single');  


max_shift_1d=ceil(2*sqrt(2)*max_shift_2d); % Maximal 1D shift that may be 
    % introduced by a 2D shift of max_shift_2d.
open_log(0);

% Find common lines.
[clstack_multi,corrstack_multi] = cryo_clmatrix_multi(pf,K,1,max_shift_1d,shift_step_1d,[],[], 5, 10);
[clstack,corrstack] = cryo_clmatrix(pf,K,1,max_shift_1d,shift_step_1d);

% Top common line between each pair of images should be the correct common
% line.
prop=comparecl( clstack_multi{1}, clstack, n_theta, thetaTol);

% Print percentage of matching common lines.
fprintf('Percentage of matching common lines = %7.4f\n',prop); % Should be very very close to 1.

correrr=norm(corrstack_multi{1}-corrstack);
fprintf('L2 difference between corrstack = %7.4f\n',correrr); % Should be very very close to 0.
