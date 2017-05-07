function [Rests,dxs,corrs]=cryo_orient_projections_gpu(projs,vol,Nrefs,trueRs,verbose,preprocess)
% CRYO_ORIENT_PROJECTION Find the orientation of a projection image.
%
% [Rests,dxs,corrs]=cryo_orient_projection(proj,vol) 
%   Given a stack of projections projs and a volume vol, estimate the
%   orientation of each projection in the stack. That is, estimate the
%   orientations in which we need to project vol to get projs.
%   the projection in the volume. Returns the estimated orientations Rest
%   and translation dxs, as well as array corrs of correlations of matches.
%   The array corrs is of dimensions size(projs,3)x2, where the i'th entry of the first column contains the correlation of the common lines between the i'th image and all reference images induced by the best matching rotation. The i'th entry of the second column contains the mean matching correlation over all tested rotations.
%
% [Rests,dxs,corrs]=cryo_orient_projections_gpu(projs,vol,Nrefs)
%   The function estimates the orientation parameters of the projections
%   projs using their common lines with Nrefs reference projections of vol.
%   By default, Nrefs is 1.5 times the dimension of vol.
%
% [Rests,dxs,corrs]=cryo_orient_projections_gpu(projs,vol,Nrefs,trueRs)
%   Provide the true rotations trueRs for debugging. Use empty array ([])
%   to ignore this parameter.
%
% [Rests,dxs,corrs]=cryo_orient_projections_gpu(projs,vol,Nrefs,trueRs,verbose)
%     0 - No screen printouts; 
%     1 - One line progress (not writteing to log); 
%     2 - Detailed progress.
%   Default is 1.
%
% [Rests,dxs]=cryo_orient_projections_gpu(projs,vol,Nrefs,trueRs,verbose,preprocess)
%   If preprocess is nonzero (default is 1) the projections and volme are
%   preprocessed. See cryo_orient_projections_auxpreprocess for details of
%   the preprocessing.
%
% To do:
%   1. Normalize rays properly (remove DC).
%   2. Filter volume and projections to the same value.
%   3. What correlation  should we expect as a function of SNR?
%
% Yoel Shkolnisky, August 2015.
% Revised: Yoel Shkolnisky, June 2016.

if ~exist('Nrefs','var') || isempty(Nrefs)
    Nrefs=-1; % Default number of references to use to orient the given 
               % projection is set below.
end

if ~exist('trueRs','var') || isempty(trueRs)
    trueRs=-1;
end

if ~exist('verbose','var')
    verbose=1;
end

if ~exist('preprocess','var')
    preprocess=1; % Default is to apply preprocessing
end

if size(projs,1)~=size(projs,2)
    error('Projections to orient must be square.');
end

szvol=size(vol);
if any(szvol-szvol(1))
    error('Volume must have all dimensions equal.');
end

if szvol(1)~=size(projs,1)
    error('Volume and projections must have matching dimensions')
end

currentsilentmode=log_silent(verbose==0);

% Preprocess projections and referece volume.
if preprocess
    log_message('Start preprocessing projections');    
    n=szvol(1);
    projs=cryo_mask(projs,1,floor(0.45*n),floor(0.05*n)); % Mask projections
    log_message('projs MD5 %s',MD5var(projs));
    log_message('Preprocessing done');
else
    log_message('Skipping preprocessing of volume and projections');
end

if Nrefs==-1
    %Nrefs=round(szvol(1)*1.5);
    Nrefs=100;
end
log_message('Using Nrefs=%d reference projections for alignment',Nrefs);

% Set the angular resolution for common lines calculations. The resolution
% L is set such that the distance between two rays that are 2*pi/L apart
% is one pixel at the outermost radius. Make L even.
% L=ceil(2*pi/atan(2/szvol(1)));
% if mod(L,2)==1 % Make n_theta even
%     L=L+1;
% end
L=360;

% Compute polar Fourier transform of the processed projections
n_r=ceil(szvol(1)/2);
log_message('Start computing polar Fourier transforms of input projections. Using n_r=%d L=%d.',n_r,L);
projs_hat=cryo_pft(projs,n_r,L,'single');
projs_hat=single(projs_hat);
log_message('Computing polar Fourier transform done');
log_message('projs_hat MD5 %s',MD5var(projs_hat));

log_message('Start normalizing Fourier transform of input projections (cryo_raynormalize)');
for k=1:size(projs,3)
    proj_hat=projs_hat(:,:,k);
% %     proj_hat=bsxfun(@times,proj_hat,H);
    proj_hat=cryo_raynormalize(proj_hat);
    projs_hat(:,:,k)=proj_hat;
end
log_message('Normalizing done');
log_message('projs_hat MD5 %s',MD5var(projs_hat));

max_shift=round(size(projs,1)*0.1);
shift_step=0.5;

[Rests,dxs,corrs]=cryo_orient_projections_gpu_worker(projs_hat,vol,Nrefs,...
    max_shift,shift_step,trueRs,verbose);

log_silent(currentsilentmode);
