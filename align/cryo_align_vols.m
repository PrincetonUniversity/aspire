function [bestR,bestdx,reflect,vol2aligned,bestcorr] = cryo_align_vols(sym,vol1,vol2,verbose,opt)
%% This function aligns vol2 according to vol1  
% Align vol2 to vol1: find the relative rotation, translation and 
% reflection between vol1 and vol2, such that vol2 is best aligned with 
% vol1.
%% input:
% sym- the symmetry type- 'Cn'\'Dn'\'T'\'O'\'I', where n is the the symmetry
%      order (for example: 'C2'). 
% vol1- 3D reference volume that vol2 should be aligned accordingly.
% vol2- 3D volume to be aligned.
% verbose- Set verbose to nonzero for verbose printouts (default is zero).
%% output: 
% bestR- the estimated rotation between vol2 and vol1, such that bestR*vol2
%        will align vol2 to vol1. 
% bestdx- size=3x1. the estimated translation between vol2 and vol1. 
% reflect- indicator for reflection. if reflect=1 then there is a
%          reflection between vol1 and vol2, else reflect=0. In order to
%          align the volumes in the case of reflect=1, you should first
%          reflect vol2 about the z axis, and then rotate by bestR.
% vol2aligned- vol2 after applyng the estimated transformation, so it is 
%              best aligned with vol1 (after optimization). 
% bestcorr- the coorelation between vol1 and vol2aligned.
%% Options:
% opt.downsample-  Downsample the volume to this size (in pixels) for
%                  faster alignment. Default is 48. Use larger value if
%                  alignment fails.
% opt.N_projs- Number of projections to use for the alignment. 
%              Defult is 30.  
% opt.true_R-  True rotation matrix between vol2 and vol1, such that 
%         vol2 =true_R*vol2. In the case of reflection, true_R should be  
%         the rotation between the volumes such that vol2=J*true_R*vol_1, 
%         where J is the reflection matrix over the z axis
%         J=diag([1,1,-1]). This input is used for debugging to calculate
%         errors.
% opt.G- Array of matrices of size 3x3xn containing the symmetry group
%        elemnts. This input is for accurate error calculation. If G is not 
%        submitted then the error will be calculated by optimization over 
%        the symmetry group 
% opt.dofscplot- set dofscplot to nonzero for FSC plot printouts (default 
%                is zero). 

% TODO (Yael, ignore this for now):
% 1. Improve printouts.

%% Check options:
defaultopt = struct('downsample',48,'N_projs',30,'G',[],'true_R',[], ...
                    'dofscplot',0);      
if ~exist('opt','var') || isempty(opt)
    opt = defaultopt;
else
    f = fieldnames(defaultopt);
    for i = 1:length(f)
        if (~isfield(opt,f{i})||(isempty(opt.(f{i}))))
            opt.(f{i}) = defaultopt.(f{i});
        end
    end
end              
%% Define variables:
N_projs = opt.N_projs; G = opt.G; true_R = opt.true_R; 
dofscplot = opt.dofscplot;
s = sym(1);
if numel(sym) > 1
        n_s = str2double(sym(2:end));
else
    n_s=0;
end
er_calc = 1;
if s == 'C' && n_s == 1, G = eye(3);
elseif isempty(G), er_calc = 0; end
refrot = 1;
if isempty(true_R), refrot = 0; end  
if ~exist('verbose','var') || isempty(verbose)
    verbose = 0;
end
currentsilentmode = log_silent(verbose == 0);
%% Validate input:
% Input volumes must be 3-dimensional, where all dimensions must be equal.
% This restriction can be remove, but then, the calculation of nr (radial
% resolution in the Fourier domain) should be adjusted accordingly. Both
% vol_1 and vol_2 must have the same dimensions.
n_1 = size(vol1);
assert(numel(n_1) == 3,'Inputs must be 3D');
assert((n_1(1) == n_1(2)),'All dimensions of input volumes must be equal');
n_2 = size(vol2);
assert(numel(n_2) == 3,'Inputs must be 3D');
assert((n_2(1) == n_1(2)) && (n_2(1) == n_1(2)),...
    'All dimensions of input volumes must be equal');
assert(n_1(1) == n_2(1),'Input volumes have different dimensions');
n = n_1(1);
n_ds = min(n,opt.downsample); % Perform aligment on down sampled volumes.  
                          % This speeds up calculation, and does not seem 
                          % to degrade accuracy
log_message('Downsampling volumes from %d to %d pixels',n,n_ds);
vol1_ds = cryo_downsample(vol1,[n_ds n_ds n_ds]);
vol2_ds = cryo_downsample(vol2,[n_ds n_ds n_ds]);
%% Aligning the volumes:
[R_est,R_est_J] = fastAlignment3D(sym,vol1_ds,vol2_ds,verbose,n_ds,N_projs,true_R,refrot,G);
vol2_aligned_ds = fastrotate3d(vol2_ds,R_est); % Rotate the original vol_2 back.
vol2_aligned_J_ds = fastrotate3d(vol2_ds,R_est_J);
vol2_aligned_J_ds = flip(vol2_aligned_J_ds,3);

estdx_ds = register_translations_3d(vol1_ds,vol2_aligned_ds);
vol2_aligned_ds = reshift_vol(vol2_aligned_ds,estdx_ds); 
estdx_J_ds = register_translations_3d(vol1_ds,vol2_aligned_J_ds);
vol2_aligned_J_ds = reshift_vol(vol2_aligned_J_ds,estdx_J_ds); 

no1 = corr(vol1_ds(:),vol2_aligned_ds(:));            
no2 = corr(vol1_ds(:),vol2_aligned_J_ds(:));          
if max(no1,no2) < 0.2 % The coorelations of the estimated rotations are 
       % smaller than 0.2, that is, no transformation was recovered. This
       % threshold was set arbitrarily.
       warning('***** Alignment failed *****');
end
%%% Do we have reflection?
reflect = 0;
corr_v = no1;
if no2 > no1
    J3 = diag([1 1 -1]); 
    corr_v = no2;
    R_est = R_est_J;
    R_est = J3*R_est*J3; 
    estdx_ds = estdx_J_ds;
    vol2_ds = flip(vol2_ds,3);  
    vol2 = flip(vol2,3);   
    reflect = 1;
    log_message('***** Reflection detected *****');
end
log_message('Correlation between downsampled aligned volumes = %7.4f',corr_v);
%% Optimization:
[bestR,~] = refine3DmatchBFGS(vol1_ds,vol2_ds,R_est,estdx_ds,0);
log_message('Done aligning downsampled volumes');
log_message('Applying estimated rotation to original volumes');
vol2aligned = fastrotate3d(vol2,bestR);
bestdx = register_translations_3d(vol1,vol2aligned);
vol2aligned = reshift_vol(vol2aligned,bestdx);    
bestcorr = corr(vol1(:),vol2aligned(:));
log_message('Estimated translations: (%7.4f,%7.4f,%7.4f)',bestdx(1),bestdx(2),bestdx(3));
log_message('Correlation between original aligned volumes = %7.4f',bestcorr);
if bestcorr < 0.5 % The coorelations of the estimated rotation are 
       % smaller than 0.5, that is, no rotation was recovered. This
       % threshold was set arbitrarily.
       warning('***** Alignment failed *****');
end
if dofscplot ~= 0
    cutoff = 0.5;
    pixA = 1;
    resA = plotFSC(vol1,vol2aligned,cutoff,pixA);
    log_message('Resolution between volumes is %5.2fA (using cutoff %4.3f)',resA,cutoff);
end
%% Error calculation:
% The difference between the estimated and reference rotation
% should be an element of the symmetry group:
if refrot ~= 0 && er_calc ~= 0
    n_g = size(G,3);
    if s == 'D' || s == 'I'
        O_g = [0 1 0; 1 0 0; 0 0 1];
        for i = 1:n_g
            G(:,:,i) = O_g*G(:,:,i)*O_g.';
        end
    end
    g_est_t = true_R.'*bestR.';
    dist = zeros(n_g,1);
    for g_idx = 1:n_g
        dist(g_idx,1) = norm(g_est_t-G(:,:,g_idx),'fro');
    end
    [~,min_idx] = min(dist);
    g_est = G(:,:,min_idx);
    err_norm = norm(true_R.' - g_est*bestR,'fro');     
    ref_true_R = true_R.';
    log_message('Reference rotation:');
    log_message('%7.4f %7.4f  %7.4f',ref_true_R(1,1),ref_true_R(1,2),ref_true_R(1,3));
    log_message('%7.4f %7.4f  %7.4f',ref_true_R(2,1),ref_true_R(2,2),ref_true_R(2,3));
    log_message('%7.4f %7.4f  %7.4f',ref_true_R(3,1),ref_true_R(3,2),ref_true_R(3,3));
    aligned_bestR = g_est*bestR;
    log_message('Estimated rotation (aligned by a symmetry element according to the reference rotation):');
    log_message('%7.4f %7.4f  %7.4f',aligned_bestR(1,1),aligned_bestR(1,2),aligned_bestR(1,3));
    log_message('%7.4f %7.4f  %7.4f',aligned_bestR(2,1),aligned_bestR(2,2),aligned_bestR(2,3));
    log_message('%7.4f %7.4f  %7.4f',aligned_bestR(3,1),aligned_bestR(3,2),aligned_bestR(3,3));    
    log_message('Estimated symmetry element:');
    log_message('%7.4f %7.4f  %7.4f',g_est(1,1),g_est(1,2),g_est(1,3));
    log_message('%7.4f %7.4f  %7.4f',g_est(2,1),g_est(2,2),g_est(2,3));
    log_message('%7.4f %7.4f  %7.4f',g_est(3,1),g_est(3,2),g_est(3,3));
    log_message('Estimation error (Frobenius norm) up to symmetry group element is %5.3e', err_norm);
    vec_ref = rotationMatrixToVector(true_R.');
    angle_ref = norm(vec_ref);
    axis_ref = vec_ref/angle_ref;    
    vec_est = rotationMatrixToVector(g_est*bestR);
    angle_est = norm(vec_est);
    axis_est = vec_est/angle_est;   
    log_message('Rotation axis:');
    log_message('\t Reference \t [%5.3f , %5.3f, %5.3f]',...
        axis_ref(1,1),axis_ref(1,2),axis_ref(1,3));
    log_message('\t Estimated \t [%5.3f , %5.3f, %5.3f]',...
        axis_est(1),axis_est(2),axis_est(3));
    log_message('Angle between axes %5.3f degrees',acosd(dot(axis_est,axis_ref)));
    log_message('In-plane rotation:');
    log_message('\t Reference \t %5.3f degrees',rad2deg(angle_ref));
    log_message('\t Estimated \t %5.3f degrees',rad2deg(angle_est));      
    log_message('\t Angle difference \t %5.3f degrees',abs(rad2deg(angle_ref)-rad2deg(angle_est))); 
end
%% Error calculation without G:
% The difference between the estimated and reference rotation
% should be an element of the symmetry group:
if refrot ~= 0 && er_calc == 0  
    G = genSymGroup(sym);   
    f = @(X)evalO(true_R.',bestR,G,X);
    opt.TolFUN = 1.0e-6;
    opt.TolX = 1.0e-6;
    opt.Display = 'off';
    X = ga(f,3,opt);
    psi = X(1); theta = X(2); phi = X(3);  
    O = Rz(phi)*Ry(theta)*Rx(psi);
    n_g = size(G,3);
    dist = zeros(1,n_g);
    for i = 1:n_g
        g = G(:,:,i);
        dist(1,i) = norm(O*true_R.' - g*O*bestR,'fro');
    end
    [err_norm,idx] = min(dist);
    g = G(:,:,idx); g_est = O.'*g*O; 
    ref_true_R = true_R.';
    log_message('Reference rotation:');
    log_message('%7.4f %7.4f  %7.4f',ref_true_R(1,1),ref_true_R(1,2),ref_true_R(1,3));
    log_message('%7.4f %7.4f  %7.4f',ref_true_R(2,1),ref_true_R(2,2),ref_true_R(2,3));
    log_message('%7.4f %7.4f  %7.4f',ref_true_R(3,1),ref_true_R(3,2),ref_true_R(3,3));
    aligned_bestR = g_est*bestR;
    log_message('Estimated rotation (aligned according to the reference by a symmetry element):');
    log_message('**** This estimation was done by optimization over the symmetry group ****');
    log_message('%7.4f %7.4f  %7.4f',aligned_bestR(1,1),aligned_bestR(1,2),aligned_bestR(1,3));
    log_message('%7.4f %7.4f  %7.4f',aligned_bestR(2,1),aligned_bestR(2,2),aligned_bestR(2,3));
    log_message('%7.4f %7.4f  %7.4f',aligned_bestR(3,1),aligned_bestR(3,2),aligned_bestR(3,3));   
    log_message('Estimated symmetry element:');
    log_message('%7.4f %7.4f  %7.4f',g_est(1,1),g_est(1,2),g_est(1,3));
    log_message('%7.4f %7.4f  %7.4f',g_est(2,1),g_est(2,2),g_est(2,3));
    log_message('%7.4f %7.4f  %7.4f',g_est(3,1),g_est(3,2),g_est(3,3));
    log_message('Estimation error (Frobenius norm) up to symmetry group element is %5.2e', err_norm);   
    vec_ref = rotationMatrixToVector(true_R.');
    angle_ref = norm(vec_ref);
    axis_ref = vec_ref/angle_ref;    
    vec_est = rotationMatrixToVector(g_est*bestR);
    angle_est = norm(vec_est);
    axis_est = vec_est/angle_est;    
    log_message('Rotation axis:');
    log_message('\t Reference \t [%5.3f , %5.3f, %5.3f]',...
        axis_ref(1,1),axis_ref(1,2),axis_ref(1,3));
    log_message('\t Estimated \t [%5.3f , %5.3f, %5.3f]',...
        axis_est(1),axis_est(2),axis_est(3));
    log_message('Angle between axes %5.2f degrees',acosd(dot(axis_est,axis_ref)));
    log_message('In-plane rotation:');
    log_message('\t Reference \t %5.3f degrees',rad2deg(angle_ref));
    log_message('\t Estimated \t %5.3f degrees',rad2deg(angle_est));      
    log_message('\t Angle difference \t %5.2f degrees',abs(rad2deg(angle_ref)-rad2deg(angle_est))); 
    log_message('**** The symmetry group of the reference volume is not pre-known and was evaluated,')
    log_message('     therefore the error calculation is not accurate ***')
end    
log_silent(currentsilentmode);
end

