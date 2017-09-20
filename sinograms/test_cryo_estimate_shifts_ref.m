% Test the function cryo_estimate_shfits_ref.
%
% The verifies the the shift_equations generated bby
% cryo_estimate_shifts_ref as well as its estimated shifts are precisely
% those estimated from cryo_clmatrix.
%
% This test code also shows how to compute rotations from quaternions, and
% how to compute a reference common lines matrix.
%
% Yoel Shkolnisky, December 2014.

n=65;
nprojs=100;
dummySNR=1;
max_shift=3;
shift_step=1;
initstate;
[projs,noisy_projs,refshifts,rots_ref]=cryo_gen_projections(n,nprojs,dummySNR,max_shift,shift_step);
masked_projs=mask_fuzzy(projs,floor(n/2)); % Applly circular mask

ref_rotations=zeros(3,3,nprojs);
 for k=1:nprojs
    ref_rotations(:,:,k)=rots_ref(:,:,k).';
 end

% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. 
n_theta=360;
n_r=ceil(n/2);
[pf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections   

%% Compte common lines matrix
% Compute reference common lines matrix
clstack_ref=clmatrix_cheat(rots_ref,n_theta);

% Search for common lines in the presence of shifts
open_log(0);
[clstack,corrstack,shift_equations,shift_equations_map,clstack_mask]=...
    cryo_clmatrix(pf,nprojs,1,ceil(2*sqrt(2)*max_shift),shift_step);

% Print the percentage of correctly detected commonlines. Shoule be very
% close to 100%.
prop=comparecl( clstack, clstack_ref, n_theta, 5 );
fprintf('Percentage of correct common lines: %f%%\n\n',prop*100);

%% Estimate rotations
S=cryo_syncmatrix_vote(clstack,n_theta);
[rotations,~,mse]=cryo_syncrotations(S,rots_ref);
fprintf('MSE of the estimated rotations: %f\n\n',mse); %Show be close to zero.

% Resgister the estimated rotations to the reference onces. 
[regrot,mse,diff,O,flag]=register_rotations(rotations,ref_rotations);
% The printed MSE should be exactly as the one two lines above.
fprintf('MSE between reference and registered estimated rotations: %f\n\n',mse);

%% Estimate shifts from equations computed by cryo_clmatrix
est_shifts=shift_equations(:,1:end-1)\shift_equations(:,end);
est_shifts=transpose(reshape(est_shifts,2,nprojs));

%% Test the function cryo_estimate_shifts_ref
% At this point we have the correct common lines matrix and the correct
% rotations.
[est_shifts_new,shift_equations_new]=cryo_estimate_shifts_ref(pf,...
    rotations,clstack,ceil(2*sqrt(2)*max_shift),shift_step);

% The following should be precisely zero.
fprintf('Difference in equations = %d\n',full(sum((shift_equations_new(:)-shift_equations(:)).^2)));
fprintf('Difference in resolving 2D shifts = %d\n',full(sum((est_shifts_new(:)-est_shifts(:)).^2)));
