function test_cryo_crop_outofcore
%
% Test the function cryo_crop_outofcore by comapring its output to
% cryo_crop. Both functions shouls have the same output.
%
% Yoel Shkolnisky, May 2016.

K=10;   % Number of projection images 
n=89;   % Each projection image is of size nxn.
max_shift_2d=0;     % Not used.
shift_step_2d=1;    % and max_shift_2d with steps of  shift_step_2d pixels.
snr=1;  % Not used. Only clean proejctions are used.

% Generate proejctions. Use only the generated clean projections "projs".
initstate;
projs=cryo_gen_projections(n,K,snr,max_shift_2d,shift_step_2d);

nc=65; % Dimensions of the cropped image.

% Crop the images in-memory
projs_cropped1=cryo_crop(projs,[nc nc],1);
projs_cropped1=single(projs_cropped1); % Convert to single precision since 
    %cryo_crop_outofcore work in single precision and we want the outputs
    %of the two functions to be compatible. 

% Crop the images from an MRC file
instackname=tempname;
outstackname=tempname;

WriteMRC(projs,1,instackname); % Save the projections to an MRC file.
cryo_crop_outofcore(instackname,outstackname,[nc nc])
projs_cropped2=ReadMRC(outstackname);

err=norm(projs_cropped1(:)-projs_cropped2(:)); % Should be zero.
fprintf('diff between two methods = %5.3e\n',err);
if err<1.0e-14
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
