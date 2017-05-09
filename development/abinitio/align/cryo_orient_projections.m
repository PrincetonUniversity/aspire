function [Rests,dxs]=cryo_orient_projections(projs,vol,Nrefs,trueRs,verbose,preprocess)
% CRYO_ORIENT_PROJECTION ind the orientation of a projection image.
%
% See cryo_orient_projections_gpu for a detailed description.
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

if size(projs,1)~=size(projs,2)
    error('Projections to orient must be square.');
end

szvol=size(vol);
if any(szvol-szvol(1))
    error('Volume must have all dimensions equal');
end

if ~exist('verbose','var')
    verbose=0;
end

currentsilentmode=log_silent(verbose==0);

% Preprocess projections and referece volume.
if ~exist('preprocess','var')
    preprocess=1; % Default is to apply preprocessing
end

if preprocess
    log_message('Preprocessing volume and projections');
    [vol,projs]=cryo_orient_projections_auxpreprocess(vol,projs);
else
    log_message('Skipping preprocessing of volume and projections');
end

szvol=size(vol); % The dimensions of vol may have changed after preprocessing.

if Nrefs==-1
    %Nrefs=round(szvol(1)*1.5);
    Nrefs=100;
end

% Generate Nrefs references projections of the given volume using random
% orientations.
log_message('Generating %d reference projections.',Nrefs);
initstate;
qrefs=qrand(Nrefs);
refprojs=cryo_project(vol,qrefs,szvol(1));
refprojs=permute(refprojs,[2 1 3]);

% Save the orientations used to generate the proejctions. These would be
% used to calculate common lines between the projection to orient and the
% reference projections.
Rrefs=zeros(3,3,Nrefs);
Rrefsvec=zeros(3,3*Nrefs);
for k=1:Nrefs
    Rrefs(:,:,k)=(q_to_rot(qrefs(:,k))).';
    Rrefsvec(:,3*(k-1)+1:3*k)=Rrefs(:,:,k);
end

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
log_message('Computing polar Fourier transforms.');
log_message('Using n_r=%d L=%d.',n_r,L);
refprojs_hat=cryo_pft(refprojs,n_r,L,'single');
projs_hat=cryo_pft(projs,n_r,L,'single');


% % if bandpass
% %     % Bandpass filter all projection, by multiplying the shift_phases the the
% %     % filter H. Then H will be applied to all projections when applying the
% %     % phases below.
% %     rk2=(0:n_r-1).';
% %     H=fuzzymask(n_r,1,floor(n_r*0.4),ceil(n_r*0.1),(n_r+1)/2).*sqrt(rk2);
% %     H=H./norm(H);    
% % else 
% %     H=1;
% % end


% Normalize polar Fourier transforms
log_message('Normalizing projections.');
for k=1:Nrefs
    pf=refprojs_hat(:,:,k);    
% %     pf=bsxfun(@times,pf,H);
    %proj(rmax:rmax+2,:)=0;
    pf=cryo_raynormalize(pf);
    refprojs_hat(:,:,k)=pf;
end

for k=1:size(projs,3)
    proj_hat=projs_hat(:,:,k);
% %     proj_hat=bsxfun(@times,proj_hat,H);
    proj_hat=cryo_raynormalize(proj_hat);
    projs_hat(:,:,k)=proj_hat;
end


% Generate candidate rotations. The rotation corresponding to the given
% projection will be searched month these rotatios.
candidate_rots=genRotationsGrid(75);
%candidate_rots(:,:,1)=Rref;
Nrots=size(candidate_rots,3);

log_message('Using %d candidate rotations.',Nrots);

% Compute the common lines between the candidate rotations and all
% reference projections. Load if possible.

log_message('Loading precomputed tables.');
%Ctbldir=fileparts(mfilename('fullpath'));
Ctbldir=tempmrcdir;
Ctblfname=fullfile(Ctbldir,'cryo_orient_projections_tables_cpu.mat');
skipprecomp=0;
if exist(Ctblfname,'file')
    tic
    precompdata=load(Ctblfname);
    t=toc;
    log_message('Loading took %5.2f seconds.',t);
    
    if isfield(precompdata,'Mkj') && ...
            isfield(precompdata,'Ckj') && isfield(precompdata,'Cjk') &&...
        isfield(precompdata,'qrefs') && isfield(precompdata,'L')
        Mkj=precompdata.Mkj;
        Ckj=precompdata.Ckj;
        Cjk=precompdata.Cjk;
        if size(Mkj,1)==Nrots && size(Mkj,2)==Nrefs && ...
                norm(qrefs-precompdata.qrefs)<1.0e-14 && ...
                L==precompdata.L
            skipprecomp=1;
        else
            log_message('Precomputed tables incompatible with input parameters.');
            log_message('\t Loaded \t Required');
            log_message('Nrots \t %d \t %d',size(Mkj,1),Nrots);
            log_message('Nrefs \t %d \t %d',size(Mkj,2),Nrefs);
            log_message('L \t \t %d \t %d',precompdata.L,L);
        end
    else
        log_message('Precomputed tables do not contain required data.');
    end
else
    log_message('Precomputed tables not found.');
end

if ~skipprecomp
    tic;
    log_message('Precomputing and saving common line tables.');
    log_message('The tables would be used in future calls to the function.');
    log_message('Patience please...');
    Ckj=(-1)*ones(Nrots,Nrefs);
    Cjk=(-1)*ones(Nrots,Nrefs);
    Mkj=zeros(Nrots,Nrefs);     % Pairs of rotations that are not "too close"
    for k=1:Nrots
        Rk=candidate_rots(:,:,k).';
        for j=1:Nrefs
            Rj=Rrefs(:,:,j).';
            if sum(Rk(:,3).*Rj(:,3)) <0.999
                %[ckj,cjk]=commonline_R(Rk,Rj,L);
                [ckj,cjk]=commonline_R2(Rk,Rj,L);
                %assert(ckj==ckj2 && cjk==cjk2);
                ckj=ckj+1; cjk=cjk+1;
                Ckj(k,j)=ckj;
                Cjk(k,j)=cjk;
                Mkj(k,j)=1;
            end
            
        end
    end
    
    t=toc;
    log_message('Precomputing tables took %5.2f seconds.',t);
        
    save(Ctblfname,'Ckj','Cjk','Mkj','qrefs', 'L');
end

% Setup shift search parameters
max_shift=round(size(projs,1)*0.1);
rmax=size(pf,1);
shift_step=0.5;
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.
log_message('Using max_shift=%d  shift_step=%d, n_shifts=%d',max_shift,shift_step,n_shifts);

rk2=(0:rmax-1).';
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(+2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
        % The shift phases are +2*pi*sqrt(-1) and not -2*pi*sqrt(-1) as
        % in cryo_estimate_shifts.m, since the phases in the latter are
        % conjugated while computing correlations, so in essesnce the
        % latter also uses  +2*pi*sqrt(-1). This function and the function 
        % cryo_estimate_shifts.m should return shifts with the same signs.            
end

% Main loop. For each candidate rotation for the givn projection, compute
% its agreement with the reference projections, and choose the best
% matching rotation.

Rests=zeros(3,3,size(projs,3));
dxs=zeros(2,size(projs,3));

t_total=tic;
parfor projidx=1:size(projs,3)
    log_message('Orienting projection %d/%d.',projidx,size(projs,3));   
    t_cpu=tic;

    cvals=zeros(Nrefs,Nrots);
    for j=1:Nrefs
        idx=find(Mkj(:,j)~=0);
        U=projs_hat(:,Ckj(idx,j),projidx);
        V=bsxfun(@times,conj(U),refprojs_hat(:,Cjk(idx,j),j));
        W=real(shift_phases.'*V);
        cvals(j,idx)=max(W);
    end
    
    scores=sum(cvals)./sum(cvals>0);
    [bestRscore,bestRidx]=max(scores);
    t_cpu=toc(t_cpu);
    log_message('\t Best correlation = %5.3f.',bestRscore);
    log_message('\t Mean correlation = %5.3f.',mean(scores));
    log_message('\t Took %5.2f seconds.',t_cpu);
    
    if trueRs~=-1
        log_message('\t Frobenius error norm of estimated rotation is %5.2e.',...
            norm(candidate_rots(:,:,bestRidx)-trueRs(:,:,projidx),'fro')/3);
    end
    
    % Now that we have the best shift, the the corresponding 2D shift of the
    % projection.
    shift_equations=zeros(Nrefs,3);
    dtheta=2*pi/L;
    
    idx=find(Mkj(bestRidx,:)==1);
    for j=idx
        pfj=bsxfun(@times,refprojs_hat(:,Cjk(bestRidx,j),j),shift_phases);
        c=real(projs_hat(:,Ckj(bestRidx,j),projidx)'*pfj);
        [~,sidx]=max(c);
        
        shift_alpha=(Ckj(bestRidx,j)-1)*dtheta;  % Angle of common ray in projection k.
        dx=-max_shift+(sidx-1)*shift_step;
        
        shift_equations(j,1)=sin(shift_alpha);
        shift_equations(j,2)=cos(shift_alpha);
        shift_equations(j,3)=dx;
    end
    dxs(:,projidx)=shift_equations(:,1:2)\shift_equations(:,3);
    Rests(:,:,projidx)=candidate_rots(:,:,bestRidx);
end
t_total=toc(t_total);
log_message('Total time for orienting %d projections is %5.2f seconds.',...
    size(projs,3),t_total);

log_silent(currentsilentmode);