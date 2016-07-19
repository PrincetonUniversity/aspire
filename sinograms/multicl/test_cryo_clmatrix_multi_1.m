% Test the function cryo_clmatrix_multi.
%
% Generate clean projections, find common lines, and compare the detected
% top common lines to the true ones.
%
% Note that if the image size is even, then we must use shift-search in
% detecting common lines, since the center is not set to half pixel (the
% center of the image) but to (0,0).
%
% Yoel Shkolnisky, May 2016.


%% Odd image, no shifts
K=10;  % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=4;     % The shifts in each projection are between -max_shift_2d  
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
shift_step_1d=0.5;    % 1D shifts steps to use in common line matching.
thetaTol=5;  % A common line is considered correct if it deviates from the 
             % true common line by up to thetaTol degrees.
snr=1;  % Not used. Only clean proejctions are used.

% Generate proejctions. Use only the generated clean projections "projs".
[projs,~,~,refq]=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

% Produce reference common lines matrix.
n_theta=360;  
[ref_clstack,~]=clmatrix_cheat_q(refq(:,1:K),n_theta);

% take Fourier transform of projections.
n_r=ceil(n/2);
pf=cryo_pft(projs,n_r,n_theta,'single');  


max_shift_1d=ceil(2*sqrt(2)*max_shift_2d); % Maximal 1D shift that may be 
    % introduced by a 2D shift of max_shift_2d.
open_log(0);

% Find common lines.
[clstack,corrstack] = cryo_clmatrix_multi(pf,K,1,max_shift_1d,shift_step_1d,[],[], 5, 10);

% Top common line between each pair of images should be the correct common
% line.
prop=comparecl( clstack{1}, ref_clstack, n_theta, thetaTol);

% Print percentage of correct common lines.
fprintf('Detection rate = %7.4f\n',prop); % Should be very very close to 1.


