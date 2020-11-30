function [estR,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA,verbose,cutoff,Rref,forcereflect)
% 
% CRYO_ALIGN_DENSITIES  Align two denisity maps
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2)
%       Align vol2 to vol1. Find the relative rotation and translation
%       between vol1 and vol2, and rotate and shift vol2 such that it is
%       best sligned with vol1. Returns the estimate rotation (Rest) and
%       tranaltion (estdx) between the volumes, as well as a copy of vol2
%       which is best aligned with vol1 (vol2aligned). The function also
%       checks if the two volumes are reflected w.r.t each other. Sets
%       reflect to 1 if reflection was detected and zero otherwise. Set 
%       verbose to nonzero for verbose printouts. 
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA)
%       Use pixel size in Angstrom pixA to compute the resolution.
%       Set to non-positive number to ignore.

% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA,verbose)
%       Set verbose to nonzero for verbose printouts. 
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,verbose,cutoff)
%       Use given FSC cutoff value. cutoff is the FSC threshold to use
%       for reporting resolutions. Default is 0.143.
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,verbose,cutoff,Rref)
%       If the true rotation between vol1 and vol2 in known (during
%       development/debugging), the function uses Rref to provide detailed
%       debugging messages. Rref is ignored if reflection is detected.
%
% Yoel Shkolnisky, December 2020.

fftw('planner','exhaustive');

if ~exist('pixA','var')
    pixA=0;
end

if ~exist('verbose','var')
    verbose=0;
end

if ~exist('cutoff','var')
    cutoff=0.143;
end

%% Validate input
% Input volumes must be 3-dimensional, where all dimensions must be equal.
% This restriction can be remove, but then, the calculation of nr (radial
% resolution in the Fourier domain) should be adjusted accordingly. Both
% vol1 and vol2 must have the same dimensions.
sz1=size(vol1);
assert(numel(sz1)==3,'Inputs must be 3D');
assert((sz1(1)==sz1(2)) && (sz1(1)==sz1(2)),...
    'All dimensions of input volumes must be equal');
sz2=size(vol2);
assert(numel(sz2)==3,'Inputs must be 3D');
assert((sz2(1)==sz1(2)) && (sz2(1)==sz1(2)),...
    'All dimensions of input volumes must be equal');
assert(sz1(1)==sz2(1),'Input volumes have different dimensions');

%% Convert reference rotation (if given) to axis/angle representatin.
% If true rotation is given (for debugging), estimate the true axis of
% rotation and the in-plane rotation angle around that axis.
%
% rotaxis_ref is the true rotation axis. gamma_ref is the in-plane rotation
% around rotaxis_ref.

refgiven=0;
if exist('Rref','var')
    if ~isempty(Rref)
        refgiven=1;
    else
        Rref=-1;
    end
else
    Rref=-1;
end

if ~exist('forcereflect','var')
    forcereflect=0;
end

%% Mask input volumes

n=size(vol1,1);
vol1masked=vol1.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));
vol2masked=vol2.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));

vol1masked=GaussFilt(vol1masked,0.3);
vol2masked=GaussFilt(vol2masked,0.3);


%% Rough alignment on downsampled volumes.
%
% Downsample the input volmes, generate a list of Nrots rotations, and find
% the rotation from this list which best aligns the two volumes. This is
% implemented by trying all rotations from the list, and for each 
% rotation looking for the best translation. The combination of
% rotation and translation which results in the best match is taken as the
% estimated relative rotation/translation between the volumes.

n_downsample=min(n,48); % Perform aligment on down sampled volumes. This 
                        % speeds up calculation, and does not seem to
                        % degrade accuracy
pixA_downsample=pixA*n/n_downsample;
vol1ds=cryo_downsample(vol1masked,[n_downsample n_downsample n_downsample]);
vol2ds=cryo_downsample(vol2masked,[n_downsample n_downsample n_downsample]);


rotations=genRotationsGrid(75);

corr0=0;
if ~forcereflect
    if verbose
        log_message('**********************************************');
        log_message('Rough alignment on downsampled masked volumes:');
        
        log_message('Volume downsampled from %d to %d',n,n_downsample);
        log_message('Using %d candidate rotations',size(rotations,3));
    end
    
    tic;
    [R0,dx0,corr0,res0,~]=bf3Dmatchaux(vol1ds,vol2ds,rotations,pixA_downsample,1,cutoff);
    t=toc;
    
    if verbose
        debugmessage(t,corr0,res0,dx0,R0,pixA_downsample,refgiven,Rref);
    end    
end

if verbose
    log_message('**************************************************************');
    log_message('Rough alignment on downsampled masked volumes with reflection:');
    
    log_message('Volume downsampled from %d to %d',n,n_downsample);
    log_message('Using %d candidate rotations',size(rotations,3));
end

tic;
[R0R,dx0R,corr0R,res0R,~]=bf3Dmatchaux(vol1ds,flip(vol2ds,3),rotations,pixA_downsample,1,cutoff);
t=toc;

assert(abs(norm(R0R)-1)<1.0e-14); % Verify that nothing bad happened
    
if verbose
    % No point in comparing to Rref if there is reflection, as Rref in this
    % case in irrelevant
    debugmessage(t,corr0R,res0R,dx0R,R0R,pixA_downsample,0,Rref);
end

reflect=0;
if corr0<corr0R
    % Reflection needed
    if verbose
        log_message('**** Reflection detected ****')
    end
    reflect=1;
    vol2ds=flip(vol2ds,3);
    R0=R0R;
    dx0=dx0R;
    
    %vol2masked=flip(vol2masked,3);
    vol2=flip(vol2,3);
end

% Now we have a rough estimate for the rotation and translation between the
% downsampled volumes. Use optimization to refine the estimates.
if verbose
    log_message('**********************************************');
    log_message('Optimized alignment of downsampled volumes:');
end

tic;
[R1,~]=refine3DmatchBFGS(vol1ds,vol2ds,R0,dx0,verbose);
t=toc;
if verbose
    vol2Rds=fastrotate3d(vol2ds,R1);
    estdx=register_translations_3d(vol1ds,vol2Rds);
    vol2Rds=reshift_vol(vol2Rds,estdx);    
    corr1=corr(vol1ds(:),vol2Rds(:));
    fsc=FSCorr(vol1ds,vol2ds);
    res=fscres(fsc,cutoff);
    pixAds=pixA*size(vol1,1)/size(vol1ds,1);
    res1=2*pixAds*numel(fsc)/res; % Resolution in Angstrom.
    debugmessage(t,corr1,res1,estdx,R1,pixA_downsample,refgiven & ~reflect,Rref);
end


% %%%%%% Debug code start
% % Check FSC between original volumes using the paramters
% % estimated thus far.
% 
% % Using R0,dx0:
% vol2Rmasked=fastrotate3d(vol2masked,R0);
% estdx=register_translations_3d(vol1masked,vol2Rmasked);
% vol2aligned=reshift_vol(vol2Rmasked,estdx);    
% plotFSC(vol1masked,vol2aligned,0.143,1);
% 
% % Using R1,dx1:
% vol2Rmasked=fastrotate3d(vol2masked,R1);
% estdx=register_translations_3d(vol1masked,vol2Rmasked);
% vol2aligned=reshift_vol(vol2Rmasked,estdx);    
% plotFSC(vol1masked,vol2aligned,0.143,1);
% %%%%%% Debug code end
if verbose
    log_message('**********************************************');
    log_message('Alignment of original (non-filtered) volumes:');
end

tic;
vol2aligned=fastrotate3d(vol2,R1);
estdx=register_translations_3d(vol1,vol2aligned);
vol2aligned=reshift_vol(vol2aligned,estdx);    
bestcorr=corr(vol1(:),vol2aligned(:));
t=toc;

estR=R1.';

fsc=FSCorr(vol1,vol2aligned);
bestRes=fscres(fsc,cutoff);
bestResA=2*pixA*numel(fsc)/bestRes; % Resolution in Angstrom.

if verbose
    debugmessage(t,bestcorr,bestResA,estdx,R1,pixA,refgiven && ~reflect,Rref);
end



function debugmessage(t,c,res,dx,R,pixA,refgiven,Rref)
% Show verbose message
log_message('Completed in %7.2f seconds',t);
log_message('Pixel size %6.3f',pixA);
log_message('Best correlation detected: %7.4f',c);
if res>0
    log_message('Best resolution detected: %5.2f',abs(res));
else
    log_message('Best resolution detected: %5.2f (dummy pixel size)',abs(res));
end
log_message('Estimated shift [%5.3f , %5.3f, %5.3f]',...
    dx(1),dx(2),dx(3));
if refgiven
    [rotaxis_ref,gamma_ref]=rot_to_axisangle(Rref);
    [rotaxis,gamma]=rot_to_axisangle(R);
    
    log_message('Rotation axis:');
    log_message('\t Estimated \t [%5.3f , %5.3f, %5.3f]',...
        rotaxis(1),rotaxis(2),rotaxis(3));
    log_message('\t Reference \t [%5.3f , %5.3f, %5.3f]',...
        rotaxis_ref(1),rotaxis_ref(2),rotaxis_ref(3));
    log_message('Angle between axes %5.2f degrees',acosd(dot(rotaxis_ref,rotaxis)));
    log_message('In-plane rotation:');
    log_message('\t Estimated \t %5.3f degrees',gamma*180/pi);
    gamma_ref_degrees=gamma_ref*180/pi;
    log_message('\t Reference \t %5.3f degrees',gamma_ref_degrees);
end
