function [Rests,dxs,corrs]=cryo_orient_projections_gpu_outofcore(projs_fname,vol,Nrefs,trueRs,verbose,preprocess)
% XXX FIX

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

projsreader=imagestackReader(projs_fname,100);
szprojs=projsreader.dim;
if szprojs(1)~=szprojs(2)
    error('Projections to orient must be square.');
end

szvol=size(vol);
if any(szvol-szvol(1))
    error('Volume must have all dimensions equal.');
end

if szvol(1)~=szprojs(1)
    error('Volume and projections must have matching dimensions')
end

currentsilentmode=log_silent(verbose==0);

processed_projs_fname=tempmrcsname;
% Preprocess projections and referece volume.
if preprocess
    log_message('Start preprocessing projections');    
    n=szvol(1);
    % Mask projections
    cryo_mask_outofcore(projs_fname,processed_projs_fname,floor(0.45*n),floor(0.05*n));
    log_message('Preprocessing done');
else
    log_message('Normalizing done');
    copyfile(projs_fname,processed_projs_fname);
    log_message('Creating copy of projections file %s',projs_fname);
end

szvol=size(vol); % The dimensions of vol may have changed after preprocessing.

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
projs_hat_fname=tempmrcsname;
log_message('Start computing polar Fourier transforms of input projections. Using n_r=%d L=%d.',n_r,L);
cryo_pft_outofcore(processed_projs_fname,projs_hat_fname,n_r,L);
log_message('Computing polar Fourier transform done');

log_message('Start normalizing Fourier transform of input projections (cryo_raynormalize)');
projs_hat_normalized_fname=tempmrcsname;
cryo_raynormalize_outofcore(projs_hat_fname,projs_hat_normalized_fname);
delete(projs_hat_fname);
log_message('Normalizing done');

max_shift=round(szprojs(1)*0.1);
shift_step=0.5;
[Rests,dxs,corrs]=cryo_orient_projections_gpu_outofcore_worker(projs_hat_normalized_fname,vol,Nrefs,...
    max_shift,shift_step,trueRs,verbose);

log_silent(currentsilentmode);

% Remove temporary files
delete(processed_projs_fname);
