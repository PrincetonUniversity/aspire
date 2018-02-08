function test_cryo_noise_estimation_outofcore
% Test the function cryo_noise_estimation_outofcore.
%
% Yoel Shkolnisky, May 2016.

K=100;   % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=0;     % Not used.
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
snr=1;  % Not used. Only clean proejctions are used.

% Generate proejctions. Use only the generated clean projections "projs".
initstate;
[~,noisy_projs]=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

% Write and read the stack to make sure we are using precisely the same
% projections.
instackname=tempname;
WriteMRC(single(noisy_projs),1,instackname); % Save the projections to an MRC file.
noisy_projs=ReadMRC(instackname);
noisy_projs=double(noisy_projs); % Although we read single precision numbers, cast to 
    % double precision since cryo_normalize_background_outofcore uses
    % double precision internally, and we want both
    % cryo_normalize_background and cryo_normalize_background_outofcore to
    % have exactly the same roundoff error.
    
% Estimate power spectrum    
psd1=cryo_noise_estimation(noisy_projs);

% Estimate noise using the images from an MRC file
psd2=cryo_noise_estimation_outofcore(instackname);

err=norm(psd1(:)-psd2(:)); % Should be zero.
fprintf('diff between two cropping methods = %5.3e\n',err); 
if err<1.0e-10
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
