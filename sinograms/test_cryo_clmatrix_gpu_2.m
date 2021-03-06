% test the function cryo_clmatrix_gpu
%
% Generate clean projections, find common lines, and compare the detected
% common lines to the true ones.
%
% Images are of even size.
%
% Yoel Shkolnisky, October 2013.
%
% Revisions:
% Y.S.  March 03, 2015    Image size set ot odd until odd/even issue in
%       gen_projections is solved.
%

K=10;  % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=2;     % The shifts in each projection are between -max_shift_2d  
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
shift_step_1d=0.5;    % 1D shifts steps to use in common line matching.
thetaTol=5;  % A common line is considered correct if it deviates from the 
             % true common line by up to thetaTol degrees.
snr=1;  % Not used. Only clean proejctions are used.

% Generate proejctions. Use only the generated clean projections "projs".
initstate;
[projs,noisy_projs,~,rots_ref]=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

% Produce reference common lines matrix.
n_theta=360;  
[ref_clstack,~]=clmatrix_cheat(rots_ref(:,:,1:K),n_theta);

% take Fourier transform of projections.
n_r=ceil(n/2);
pf=cryo_pft(projs,n_r,n_theta,'single');  


max_shift_1d=ceil(2*sqrt(2)*max_shift_2d); % Maximal 1D shift that may be 
    % introduced by a 2D shift of max_shift_2d.
max_shift_1d=max(max_shift_1d,2); % For an even image we must test for 
    % shifts of at least 2 pixels since PFT does not use a half pixel as a 
    % center. If you use max_shift_1d=0, the detection rate drops to 91% 
    % (for clean images)
open_log(0);

% Find common lines.
clstack = cryo_clmatrix_gpu(pf,K,1,max_shift_1d,shift_step_1d);
prop=comparecl( clstack, ref_clstack, n_theta, thetaTol);

% Print percentage of correct common lines.
fprintf('Detection rate = %7.4f\n',prop); % Should be very close to 1.
