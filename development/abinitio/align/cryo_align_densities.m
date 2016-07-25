function [estR,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA,verbose,Rref,forcereflect,Nprojs)
% CRYO_ALIGN_DENSITIES  Align two denisity maps
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2)
%       Align vol2 to vol1. Find the relative rotation and translation
%       between vol1 and vol2, and rotate and shift vol2 such that it is
%       best aligned with vol1. Returns the estimated rotation (Rest) and
%       tranaltion (estdx) between the volumes, as well as a copy of vol2
%       which is best aligned with vol1 (vol2aligned). The function also
%       checks if the two volumes are reflected w.r.t each other. Sets
%       reflect to 1 if reflection was detected (zero otherwise). 
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA)
%       Use pixel size in Angstrom pixA to compute the resolution.
%       Set to non-positive number to ignore.
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA)
%       Set verbose to nonzero for verbose printouts. 
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,pixA,verbose)
%       Set verbose to nonzero for verbose printouts (default is zero).
%
% [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(vol1,vol2,verbose,Rref)
%       Use the true rotation between vol1 and vol2 specified by Rref to
%       provide detailed debugging messages. Rref is ignored if reflection
%       is detected.
%
% Yoel Shkolnisky, June 2016.


if ~exist('pixA','var')
    pixA=0;
end

if ~exist('verbose','var')
    verbose=0;
end

if ~exist('Nprojs','var')
    Nprojs=100;  % Number of projections to use for alignment.
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

timing=tic;

%% Mask input volumes

n=size(vol1,1);
vol1masked=vol1.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));
vol2masked=vol2.*fuzzymask(n,3,floor(0.45*n),floor(0.05*n));

vol1masked=GaussFilt(vol1masked,0.3);
vol2masked=GaussFilt(vol2masked,0.3);


%% Alignment 

refq=qrand(Nprojs); % Quaternions used to project volume 2.

% Convert quaternions to rotations
trueRs=zeros(3,3,Nprojs); 
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(refq(:,k))).';
end

% Generate projections of volume 2.
projs2=cryo_project(vol2masked,refq);
projs2=permute(projs2,[2,1,3]);

% Estimate rotations of the projections of volume 2.
log_message('Aligning volumes.')
[Rests,dxests]=cryo_orient_projections_gpu(projs2,vol1masked,Nprojs,[],verbose,0);

% Assess quality of the alignment. Use Rests and trueRs to estimate the
% matrix aligning the two volumes. The close this matrix to an orthogonal
% matrix, the better the alignment. 
%
% There are two possibilties:
% 1. There is no reflection between vol1 and vol2, in which case they are
% related by rotation only, that is R_{i} = O \tilde{R}_{i}, where R_{i}
% and tilde{R}_{i} are the rotations are the rotations of the projections
% in the coordinates systems of vol2 and vol1, respectively (corresponding
% to trueRs and Rests, respectively).
% 2. There is reflection, in which case R_{i} = O J \tilde{R}_{i} J.
%
% We estimate O in both cases and choose the one for which the resulting O
% is more orthogonal.
R1=zeros(3);
Omat=zeros(3);
J3=diag([1 1 -1]);
for k=1:Nprojs
    R1=R1+trueRs(:,:,k)*Rests(:,:,k).';
    Omat=Omat+trueRs(:,:,k)*(J3*Rests(:,:,k)*J3).';
end
R1=R1./Nprojs;
Omat=Omat./Nprojs;

s1=svd(R1); s2=svd(Omat);
no1=max(s1)/min(s1); %Non orthogonality (actually, the condition number)
no2=max(s2)/min(s2);

if verbose
    log_message('Singular values of aligning matrix:');
    log_message('\t without reflection (%7.4f,%7.4f,%7.4f); condition number = %7.4f',s1(1),s1(2),s1(3),no1);
    log_message('\t with    reflection (%7.4f,%7.4f,%7.4f); condition number = %7.4f',s2(1),s2(2),s2(3),no2);
end


log_message('Refining alignment.');
Rests=cryo_refine_orientations(projs2,vol1masked,Rests,dxests,1,-1);


% There are two possibilties:
% 1. There is no reflection between vol1 and vol2, in which case they are
% related by rotation only, that is R_{i} = O \tilde{R}_{i}, where R_{i}
% and tilde{R}_{i} are the rotations are the rotations of the projections
% in the coordinates systems of vol2 and vol1, respectively (corresponding
% to trueRs and Rests, respectively).
% 2. There is reflection, in which case R_{i} = O J \tilde{R}_{i} J.
%
% We estimate O in both cases and choose the one for which the resulting O
% is more orthogonal.
R1=zeros(3);
Omat=zeros(3);
J3=diag([1 1 -1]);
for k=1:Nprojs
    R1=R1+trueRs(:,:,k)*Rests(:,:,k).';
    Omat=Omat+trueRs(:,:,k)*(J3*Rests(:,:,k)*J3).';
end
R1=R1./Nprojs;
Omat=Omat./Nprojs;

s1=svd(R1); s2=svd(Omat);
no1=max(s1)/min(s1); %Non orthogonality (actually, the condition number)
no2=max(s2)/min(s2);

if verbose
    log_message('Singular values of aligning matrix:');
    log_message('\t without reflection (%7.4f,%7.4f,%7.4f); condition number = %7.4f',s1(1),s1(2),s1(3),no1);
    log_message('\t with    reflection (%7.4f,%7.4f,%7.4f); condition number = %7.4f',s2(1),s2(2),s2(3),no2);
end

if min(no1,no2)>1.2 % The condition number of the estimated rotation is 
       % larger than 1.2, that is, no rotation was recovered. This
       % threshold was set arbitrarily.
       warning('Alignment failed.');
end

reflect=0; % Do we have reflection?
R=R1;
if no2<no1 || forcereflect
    reflect=1;
    R=Omat;
    
    if verbose
        log_message('**** Reflection detected ****');
    end
end
    
[U,~,V]=svd(R); % Project R to the nearest rotation.
Omat=U*V.';
assert(det(Omat)>0);
Omat=Omat([2 1 3],[2 1 3]);


if verbose
    log_message('Estimated rotation:');
    log_message('%7.4f %7.4f  %7.4f',Omat(1,1),Omat(1,2),Omat(1,3));
    log_message('%7.4f %7.4f  %7.4f',Omat(2,1),Omat(2,2),Omat(2,3));
    log_message('%7.4f %7.4f  %7.4f',Omat(3,1),Omat(3,2),Omat(3,3));
    
    if refgiven
        tmp=Rref;
        
        if reflect
            %tmp=J3*tmp*J3;
            % Seems it should be tmp=J2*tmp*J3, but I did not look into
            % that.
        end
        log_message('Reference rotation:');
        log_message('%7.4f %7.4f  %7.4f',tmp(1,1),tmp(1,2),tmp(1,3));
        log_message('%7.4f %7.4f  %7.4f',tmp(2,1),tmp(2,2),tmp(2,3));
        log_message('%7.4f %7.4f  %7.4f',tmp(3,1),tmp(3,2),tmp(3,3));
        if reflect
            %log_message('Reference rotation was J-conjugated');
        end
        
        log_message('Estimation error (Frobenius norm) = %5.3e', norm(Omat-tmp,'fro'));
    end
end


vol2maskedaligned=fastrotate3d(vol2masked,Omat.'); % Rotate the masked vol2 back.
vol2aligned=fastrotate3d(vol2,Omat.'); % Rotate the original vol2 back.

if reflect
    vol2maskedaligned=flip(vol2maskedaligned,3);
    vol2aligned=flip(vol2aligned,3);
end

estdx=register_translations_3d(vol1masked,vol2maskedaligned);
if verbose
    log_message('Estimated translations: (%7.4f,%7.4f,%7.4f)',estdx(1),estdx(2),estdx(3));
end

vol2maskedaligned=reshift_vol(vol2maskedaligned,estdx);
vol2aligned=reshift_vol(vol2aligned,estdx);    

timing=toc(timing);

estR=Omat;

% Compute correlations
c_masked=corr(vol1masked(:),vol2maskedaligned(:)); % Masked volumes
c_orig=corr(vol1(:),vol2aligned(:)); % Original volumes
fsc=FSCorr(vol1,vol2aligned);
res=fscres(fsc,0.134);
resA=2*pixA*numel(fsc)/res; % Resolution in Angstrom.

log_message('Completed in %7.2f seconds',timing);
log_message('Pixel size %6.3f',pixA);
log_message('Correlation between masked aligned volumes = %7.4f',c_masked);
log_message('Correlation between original aligned volumes = %7.4f',c_orig);

if res>0
    log_message('Resolution of aligned volumes = %5.2f',abs(resA));
else
    log_message('Resolution of aligned volumes = %5.2f (dummy pixel size)',abs(resA));
end

if refgiven
    [rotaxis_ref,gamma_ref]=rot_to_axisangle(Rref);
    [rotaxis,gamma]=rot_to_axisangle(Omat.');
    
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
