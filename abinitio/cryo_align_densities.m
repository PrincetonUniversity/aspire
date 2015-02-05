function [estR,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA,verbose,Rref)
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

% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA)
%       Set verbose to nonzero for verbose printouts. 

% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,verbose,Rref)
%       If the true rotation between vol1 and vol2 in known (during
%       development/debugging), the function uses Rref to provide detailed
%       debugging messages. Rref is ignored if reflection is detected.
%
% Yoel Shkolnisky, January 2015.

if ~exist('pixA','var')
    pixA=0;
end

if ~exist('verbose','var')
    verbose=0;
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
    refgiven=1;
else
    Rref=-1;
end

%% Mask input volumes

n=size(vol1,1);
vol1=vol1.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));
vol2=vol2.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));


%% Rough alignment on downsampled volumes.
%
% Downsample the input volmes, generate a list of Nrots rotations, and find
% the rotation from this list which best aligns the two volumes. This is
% implemented by trying all rotations from the list, and for each 
% rotation looking for the best translation. The combination of
% rotation and translation which results in the best match is taken as the
% estimated relative rotation/translation between the volumes.

n_downsample=33;
pixA_downsample=pixA*n/n_downsample;
vol1ds=Downsample(vol1,[n_downsample n_downsample n_downsample]);
vol2ds=Downsample(vol2,[n_downsample n_downsample n_downsample]);


Nrots=5000;
qrots=qrand(Nrots);
rotations=zeros(3,3,Nrots);
for jj=1:Nrots
    rotations(:,:,jj)=q_to_rot(qrots(:,jj));
end

if verbose
    fprintf('**********************************************\n');
    fprintf('Rough alignment on downsampled masked volumes:\n');
end

tic;
[R0,dx0,corr0,res0,~]=bf3Dmatchaux(vol1ds,vol2ds,rotations,pixA_downsample,1);
t=toc;

if verbose
    debugmessage(t,corr0,res0,dx0,R0,refgiven,Rref);
end


if verbose
    fprintf('**************************************************************\n');
    fprintf('Rough alignment on downsampled masked volumes with reflection:\n');
end

tic;
[R0R,dx0R,corr0R,res0R,~]=bf3Dmatchaux(vol1ds,flipdim(vol2ds,3),rotations,pixA_downsample,1);
t=toc;

if verbose
    % No point in comparing to Rref if there is reflection, as Rref in this
    % case in irrelevant
    debugmessage(t,corr0R,res0R,dx0R,R0R,0,Rref);
end

reflect=0;
if corr0<corr0R
    % Reflection needed
    fprintf('**** Reflection detected ****.\n')
    reflect=1;
    vol2ds=flipdim(vol2ds,3);
    R0=R0R;
end

%% Refine search on downsampled volumes
%
% Take the estimated relative rotation/translation found above and refine
% it. This is implemented by selecting several random rotations around the
% estimated rotation, and for each such rotation looking for the best
% translation. The best match is taken as the new estimate of the relative
% rotation/translation.

if verbose
    fprintf('************************************************\n');
    fprintf('Refined alignment on downsampled masked volumes:\n');
end
   
newrots=genNearRotations(R0,10,100,10,31);

tic;
[R1,dx1,corr1,res1,~]=bf3Dmatchaux(vol1ds,vol2ds,newrots,pixA_downsample,1);
t=toc;


if verbose
    debugmessage(t,corr1,res1,dx1,R1,refgiven && ~reflect,Rref);
end

%% Refine once more on the original volume
%
% Use the prevoius estimated rotation/translation as a starting point for
% aligning the two original (not downsampled) volumes. This is done as
% above by using several random rotation near the previously estimated
% rotation.

if verbose
    fprintf('*********************************************\n');
    fprintf('Refined alignment on original masked volumes:\n');
end

newrots=genNearRotations(R1,5,50,5,31);

if reflect
    vol2=flipdim(vol2,3);
end

tic;
[bestR,bestdx,bestcorr,bestRes,~]=bf3Dmatchaux(vol1,vol2,newrots,pixA,1);
t=toc;

if verbose
    debugmessage(t,bestcorr,bestRes,bestdx,bestR,refgiven && ~reflect,Rref);
end

vol2aligned=fastrotate3d(vol2,bestR);
vol2aligned=reshift_vol(vol2aligned,bestdx);    
estR=bestR.';
estdx=bestdx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotations=genNearRotations(R,axisdeviation,Naxes,gammadeviation,Ngamma)

% Create rotations near the estimated rotation. Find a plane perpendicular
% to the estimated rotation axis, find points close to the rotation axis on
% this plane, and project them back to the spehere. These are the axes of
% the new rotations to try.

[rotaxis,gamma]=rot_to_axisangle(R);

[Q,~]=qr([rotaxis rand(3,2)]);
assert(norm(Q*Q.'-eye(3))<1.0e-14,'Q not orthogonal'); 

if abs(det(Q)-1)>0.1 % Vector are left handed so flip
    Q(:,[3,2])=Q(:,[2,3]);
end

% Basis for the plane perpendicular to the projection direction rotaxis.
u_x=Q(:,2);
u_y=Q(:,3);

%Naxes=50;  % Generates Naxes new rotation axes
alpha=(rand(2,Naxes)-1/2)*2; % Random coordinates in [-1,1] of points on 
                             % the plane perpendicular to rotaxis.
alpha(:,1)=[0 0].'; % Make the current best guess one of the points
% 0.1 corresponds roughly to a deviation of +/- 5 degrees from the rotation
% axis
newaxes=repmat(rotaxis,1,Naxes)...
    +sin(axisdeviation).*repmat(u_x,1,Naxes).*repmat(alpha(1,:),3,1)...
    +sin(axisdeviation).*repmat(u_y,1,Naxes).*repmat(alpha(2,:),3,1);
newaxes=newaxes./repmat(sqrt(sum(newaxes.^2,1)),3,1); % Project the axes back to the sphere.

% Generate new in-planes rotations that deviate by +/- degrees from the
% estimated in-plane rotation
gamma_degress=gamma/pi*180;
if mod(Ngamma,2)==0
    Ngamma=Ngamma+1; % Make sure to have an odd number of points so that 
                     % the current best guess is one of the points.
end
newgammas=(linspace(-gammadeviation,gammadeviation,Ngamma)+gamma_degress)/180*pi;

rotations=zeros(3,3,Naxes*numel(newgammas));
idx=1;
for kk=1:Naxes
    for jj=1:numel(newgammas)
        rotations(:,:,idx)=axisangle_to_rot(newaxes(:,kk),newgammas(jj));
        idx=idx+1;
    end
end


function debugmessage(t,c,res,dx,R,refgiven,Rref)
% Show verbose message
fprintf('Completed in %7.2f seconds\n',t);
fprintf('Best correlation detected: %7.4f\n',c);
fprintf('Best resolution detected: %5.2f',abs(res));
if res<0
    fprintf(' (dummy pixel size)');
end
fprintf('\n');
fprintf('Estimated shift [%5.3f , %5.3f, %5.3f]\n',...
    dx(1),dx(2),dx(3));
if refgiven
    [rotaxis_ref,gamma_ref]=rot_to_axisangle(Rref);
    [rotaxis,gamma]=rot_to_axisangle(R);
    
    fprintf('Rotation axis:\n');
    fprintf('\t Estimated \t [%5.3f , %5.3f, %5.3f]\n',...
        rotaxis(1),rotaxis(2),rotaxis(3));
    fprintf('\t Reference \t [%5.3f , %5.3f, %5.3f]\n',...
        rotaxis_ref(1),rotaxis_ref(2),rotaxis_ref(3));
    fprintf('Angle between axes %5.2f degrees\n',acosd(dot(rotaxis_ref,rotaxis)));
    fprintf('In-plane rotation:\n');
    fprintf('\t Estimated \t %5.3f degrees\n',gamma*180/pi);
    gamma_ref_degrees=gamma_ref*180/pi;
    fprintf('\t Reference \t %5.3f degrees\n',gamma_ref_degrees);
end