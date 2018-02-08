function test_cryo_mask_outofcore%
% Test the function cryo_mask_outofcore by comapring its output to
% cryo_mask. Both functions shouls have the same output.
%
% Yoel Shkolnisky, April 2017.

K=10;   % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=0;     % Not used.
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
snr=1;  % Not used. Only clean proejctions are used.

% Generate proejctions. Use only the generated clean projections "projs".
initstate;
projs=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

r=floor(0.45*n);
risetime=floor(0.05*n);
% Crop the images in-memory
projs=single(projs); % Convert to single precision since 
    %cryo_mask_outofcore work in single precision and we want the outputs
    %of the two functions to be compatible. 
projs_masked1=cryo_mask(projs,1,r,risetime);

% Crop the images from an MRC file
instackname=tempname;
outstackname=tempname;

WriteMRC(projs,1,instackname); % Save the projections to an MRC file.
cryo_mask_outofcore(instackname,outstackname,r,risetime)
projs_masked2=ReadMRC(outstackname);

err=norm(projs_masked1(:)-projs_masked2(:)); % Should be zero.
fprintf('diff between two methods = %5.3e\n',err);
if err<1.0e-14
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
