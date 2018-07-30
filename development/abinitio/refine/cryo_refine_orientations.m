function [rots_refined,shifts_refined,errvec]=cryo_refine_orientations(projs,isFT,vol,rots,shifts,verbose,Nrefs,true_Rs,true_shifts)

% if FT is 0 if we get real space projections and 1 if we get polar-Fourier
% transformed projections

%(proj_hat,R,refprojs_hat,Rrefs,L,estdx)

if ~exist('Nrefs','var') || isempty(Nrefs)
    Nrefs=-1; % Default number of references to use to orient the given 
               % projection is set below.
end

if ~exist('true_Rs','var') || isempty(true_Rs)
    true_Rs=-1;
end

if ~exist('true_shifts','var') || isempty(true_shifts)
    true_shifts=-1;
end

szvol=size(vol);
if any(szvol-szvol(1))
    error('Volume must have all dimensions equal');
end

if ~isFT
    if size(projs,1)~=size(projs,2)
        error('Projections to orient must be square.');
    end
    
    if size(projs,1)~=szvol(1)
        error('Projections (projs) and volume (vol) must has same dimensions.');
    end
end

if ~exist('verbose','var')
    verbose=1;
end

Nprojs=size(projs,3);
if size(rots,3)~=Nprojs
    error('Number of rotations (rots) must be equal to the number of projections.');
end

if size(shifts,2)~=Nprojs
     error('Number of shifts must be equal to the number of projections.');
end

currentsilentmode=log_silent(verbose==0);

if Nrefs==-1
    %Nrefs=round(szvol(1)*1.5);
    Nrefs=100;
end

if isFT
    % The Fourier transformsed projecitons are assumed to be normalized
    % using cryo_raynormalize. This is not enforced here to save time with
    % large datasets.
    projs_hat=projs;
else
    % Set the angular resolution for common lines calculations. The resolution
    % L is set such that the distance between two rays that are 2*pi/L apart
    % is one pixel at the outermost radius. Make L even.
    % L=ceil(2*pi/atan(2/szvol(1)));
    % if mod(L,2)==1 % Make n_theta even
    %     L=L+1;
    % end
    L=360;
    % Compute polar Fourier transform of the projecitons.
    n_r=ceil(szvol(1)/2);
    projs_hat=cryo_pft(projs,n_r,L,'single'); % XXX No reason to recompute that. Just get it as parameter.
    projs_hat=single(projs_hat);
    projs_hat=cryo_raynormalize(projs_hat);
end
log_message('projs_hat MD5 %s',MD5var(projs_hat));

% Generate reference projections
rots_ref=rand_rots(Nrefs);  % Generate Nprojs projections to orient.
projs_ref=cryo_project(vol,rots_ref);
projs_ref=permute(projs_ref,[2,1,3]);
Rrefs = permute(rots_ref, [2 1 3]);

% Compute polar Fourier transform of reference projections.
n_r=size(projs_hat,1);
L=size(projs_hat,2);
projs_ref_hat=cryo_pft(projs_ref,n_r,L,'single');
log_message('projs_ref_hat MD5 %s',MD5var(projs_ref_hat));
projs_ref_hat=cryo_raynormalize(projs_ref_hat);
log_message('projs_ref_hat MD5 %s',MD5var(projs_ref_hat));

rots_refined=zeros(size(rots));
shifts_refined=zeros(size(shifts));
errvec=zeros(Nprojs,4);

printProgressBarHeader;
parfor k=1:Nprojs
    progressTic(k,Nprojs);
    dx=shifts(:,k); 
    R=rots(:,:,k);
    [rots_refined(:,:,k),shifts_refined(:,k),~]=optimize_orientation(projs_hat(:,:,k),R,projs_ref_hat,Rrefs,L,dx);

%     if (numel(true_Rs)>1) || (numel(true_shifts)>1)            
%         fprintf('k=%d/%d\n',k,Nprojs);
%     end
    

    % Measure errors relative to true orientation parameters.
    rot_L2_error_before_refinement=0;
    rot_L2_error_after_refinement=0;
    shifts_L2_error_before_refinement=0;
    shifts_L2_error_after_refinement=0;
    
    if numel(true_Rs)>1
        rot_L2_error_before_refinement=norm(R-true_Rs(:,:,k),'fro')/norm(true_Rs(:,:,k),'fro');
        rot_L2_error_after_refinement=norm(rots_refined(:,:,k)-true_Rs(:,:,k),'fro')/norm(true_Rs(:,:,k),'fro');
        %fprintf('\t Rotation error: before=%e \t after=%e\n',rot_L2_error_before_refinement,rot_L2_error_after_refinement);
    end
    
    if numel(true_shifts)>1
        shifts_L2_error_before_refinement=norm(shifts(:,k)-true_shifts(:,k));
        shifts_L2_error_after_refinement=norm(shifts_refined(:,k)-true_shifts(:,k));
        %fprintf('\t Shift error:    before=%e \t after=%e\n',shifts_L2_error_before_refinement,shifts_L2_error_after_refinement);
    end

    errvec(k,:)=[rot_L2_error_before_refinement ...
        rot_L2_error_after_refinement ...
        shifts_L2_error_before_refinement ...
        shifts_L2_error_after_refinement];
end

log_silent(currentsilentmode);
