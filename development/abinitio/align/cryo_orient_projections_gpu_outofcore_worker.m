function [Rests,dxs,corrs]=cryo_orient_projections_gpu_outofcore_worker(projs_hat_normalized_fname,vol,Nrefs,...
    max_shift,shift_step,trueRs,verbose)

szvol=size(vol);
n=szvol(1);
vol=cryo_mask(vol,0,floor(0.45*n),floor(0.05*n)); % Mask volume
log_message('vol MD5 %s',MD5var(vol));

projshat_reader=imagestackReaderComplex(projs_hat_normalized_fname);
szprojs=projshat_reader.dim;
n_r=szprojs(1);
L=szprojs(2);
n_projs=szprojs(3);


% Generate Nrefs references projections of the given volume using random
% orientations.
log_message('Start generating %d reference projections of size %dx%d.',Nrefs,szvol(1),szvol(1));
qrefs=qrand(Nrefs);
refprojs=cryo_project(vol,qrefs,szvol(1));
refprojs=permute(refprojs,[2 1 3]);
log_message('Generating reference projections done');
log_message('qrefs MD5 %s',MD5var(qrefs));
log_message('refprojs MD5 %s',MD5var(refprojs));

% Save the orientations used to generate the proejctions. These would be
% used to calculate common lines between the projection to orient and the
% reference projections.
Rrefs=zeros(3,3,Nrefs);
Rrefsvec=zeros(3,3*Nrefs);
for k=1:Nrefs
    Rrefs(:,:,k)=(q_to_rot(qrefs(:,k))).';
    Rrefsvec(:,3*(k-1)+1:3*k)=Rrefs(:,:,k);
end
 
% Compute polar Fourier transform of reference projecitons.
log_message('Start computing polar Fourier transforms of reference projections. Using n_r=%d L=%d.',n_r,L);
refprojs_hat=cryo_pft(refprojs,n_r,L,'single');
log_message('Computing polar Fourier transform done');
log_message('refprojs_hat MD5 %s',MD5var(refprojs_hat));

% Normalize polar Fourier transforms
log_message('Start normalizing Fourier transform of reference projections (cryo_raynormalize)');
for k=1:Nrefs
    pf=refprojs_hat(:,:,k);    
% %     pf=bsxfun(@times,pf,H);
    %proj(rmax:rmax+2,:)=0;
    pf=cryo_raynormalize(pf);
    refprojs_hat(:,:,k)=pf;
end
log_message('Normalizing done');
log_message('refprojs_hat MD5 %s',MD5var(refprojs_hat));

% Generate candidate rotations. The rotation corresponding to the given
% projection will be searched month these rotatios.
candidate_rots=genRotationsGrid(75);
log_message('candidate_rots MD5 %s',MD5var(candidate_rots));
%candidate_rots=genRotationsGrid(50);
%candidate_rots(:,:,1)=Rref;
Nrots=size(candidate_rots,3);

log_message('Using %d candidate rotations for alignment.',Nrots);

% Compute the common lines between the candidate rotations and all
% reference projections. Load if possible.

log_message('Loading precomputed tables.');
%Ctbldir=fileparts(mfilename('fullpath'));
Ctbldir=tempmrcdir;
Ctblfname=fullfile(Ctbldir,'cryo_orient_projections_tables_gpu.mat');
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
    log_message('Precomputing and saving common line tables.');
    log_message('The tables would be used in future calls to the function.');
    log_message('Patience please...');

%     % Slower implementation of the code below.
%     % To compare between the two code change in the slow code:
%     %     Ckj  to  C2kj
%     %     Cjk  to  C2jk
%     %     Mkj  to  M2kj
%     tic;
%     Ckj=(-1)*ones(Nrots,Nrefs);
%     Cjk=(-1)*ones(Nrots,Nrefs);
%     Mkj=zeros(Nrots,Nrefs);     % Pairs of rotations that are not "too close"
%     for k=1:Nrots
%         Rk=candidate_rots(:,:,k).';
%         for j=1:Nrefs
%             Rj=Rrefs(:,:,j).';
%             if sum(Rk(:,3).*Rj(:,3)) <0.999
%                 [ckj,cjk]=commonline_R2(Rk,Rj,L);
%                 ckj=ckj+1; cjk=cjk+1;
%                 Ckj(k,j)=ckj;
%                 Cjk(k,j)=cjk;
%                 Mkj(k,j)=1;
%             end
%             
%         end
%     end
%     toc
%     % END slower implementation
    tic;
    
    candidate_rots_vec=zeros(3,3*Nrots);
    candidate_rots_t_vec=zeros(3*Nrots,3);
    for k=1:Nrots
        candidate_rots_vec(:,3*(k-1)+1:3*k)=candidate_rots(:,:,k);
        candidate_rots_t_vec(3*(k-1)+1:3*k,:)=(candidate_rots(:,:,k)).';
    end
    
    Ckj=(-1)*ones(Nrots,Nrefs);
    Cjk=(-1)*ones(Nrots,Nrefs);
    Mkj=zeros(Nrots,Nrefs);     % Pairs of rotations that are not "too close"
       
    for j=1:Nrefs        
        Rj=Rrefs(:,:,j);       
        [Mkj(:,j),Ckj(:,j),Cjk(:,j)]=commonline_R_vec(candidate_rots,Rj,L,0.999);     
    end
    
    % Convert to single to save space in the tables
    Ckj=single(Ckj);
    Cjk=single(Cjk);
    Mkj=single(Mkj);
    
    t=toc;
    log_message('Precomputing tables took %5.2f seconds.',t);
        
    save(Ctblfname,'Ckj','Cjk','Mkj','qrefs', 'L');
    system(sprintf('chmod a+rwx %s',Ctblfname));
end

% Setup shift search parameters
rmax=n_r;
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

Rests=zeros(3,3,n_projs);
dxs=zeros(2,n_projs);


gCkj=gpuArray(single(Ckj));
gCjk=gpuArray(single(Cjk));
grefprojs_hat=gpuArray(single(refprojs_hat));
gshift_phases=gpuArray(single(shift_phases));


% Precompute all pointers for extracting commong lines.
idx_cell=cell(Nrefs,1);
gCkj_cell=cell(Nrefs,1);    
gCjk_cell=cell(Nrefs,1);    
grefprojs_cell=cell(Nrefs,1);
for j=1:Nrefs
        idx=find(Mkj(:,j)~=0);
        idx_cell{j}=idx;
        gidx=gpuArray(single(idx));
        gCkj_cell{j}=gCkj(gidx,j);
        gCjk_cell{j}=gCjk(gidx,j);
        grefprojs_cell{j}=grefprojs_hat(:,gCjk_cell{j},j);
end
gshift_phases=gshift_phases.';

corrs=zeros(n_projs,2); % Statistics on common-lines matching.

progress_msg=[]; % String use to print progress message.

t_total=tic;
for projidx=1:szprojs(3)
    if verbose==2
    	log_message('Orienting projection %d/%d.',projidx,n_projs);   
    end
    t_gpu=tic;
    currproj=projshat_reader.getImage(projidx);
    gproj_hat=gpuArray(single(currproj));   
    cvals=zeros(Nrefs,Nrots);
    for j=1:Nrefs
%         idx=find(Mkj(:,j)~=0);
%         gidx=gpuArray(single(idx));
%         gU=gproj_hat(:,gCkj(gidx,j));
%         gV=bsxfun(@times,conj(gU),grefprojs_hat(:,gCjk(gidx,j),j));
%         gW=real(gshift_phases.'*gV);
%         cvals(j,idx)=gather(max(gW));
        
%         % Optimized version of the above code.
%         gU=gproj_hat(:,gCkj_cell{j});
%         gV=bsxfun(@times,conj(gU),grefprojs_hat(:,gCjk_cell{j},j));
%         gW=real(gshift_phases.'*gV);
%         cvals(j,idx_cell{j})=gather(max(gW));
        
        % One more optimization:
        gU=gproj_hat(:,gCkj_cell{j});
        gV=bsxfun(@times,conj(gU),grefprojs_cell{j});
        gW=real(gshift_phases*gV);
        cvals(j,idx_cell{j})=gather(max(gW));

    end
    
    scores=sum(cvals)./sum(cvals>0);
    [bestRscore,bestRidx]=max(scores);
    meanRscore=mean(scores);
    
    t_gpu=toc(t_gpu);
    if verbose==1
        if mod(projidx,50)==0
            bs=char(repmat(8,1,numel(progress_msg)));
            fprintf('%s',bs);
            progress_msg=sprintf(...
                'Orienting projection %d/%d (corr:best=%7.4f, mean=%7.4f, b/m=%7.2f)  t=%5.2f secs',...
                projidx,n_projs,bestRscore,meanRscore,bestRscore/meanRscore,t_gpu);
            fprintf('%s',progress_msg);
        end
        
    elseif verbose==2
        log_message('\t Best correlation = %5.3f.',bestRscore);
        log_message('\t Mean correlation = %5.3f.',meanRscore);
        log_message('\t Took %5.2f seconds.',t_gpu);
    end
    
    corrs(projidx,1)=bestRscore;
    corrs(projidx,2)=meanRscore;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % Non GPU code - uncomment if no GPU exist.
    %     tic;
    %     cvals=zeros(Nrefs,Nrots);
    %     for j=1:Nrefs
    %         idx=find(Mkj(:,j)~=0);
    %         U=projs_hat(:,Ckj(idx,j),projidx);
    %         V=bsxfun(@times,conj(U),refprojs_hat(:,Cjk(idx,j),j));
    %         W=real(shift_phases.'*V);
    %         cvals(j,idx)=max(W);
    %     end
    %     % To run comparison against refernce and GPU code rename in the
    %     % following line
    %     %  bestRscore to bestRscore_cpu
    %     %  bestRidx   to bestRidx_cpu
    %    scores=sum(cvals)./sum(cvals>0);
    %    [bestRscore,bestRidx]=max(scores);
    %     t_cpu=toc;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Rerefence code
    % tic
    % parfor k=1:Nrots
    %     cvals=zeros(Nrefs,1);
    %
    %     idx=find(Mkj(k,:)==1);
    %     for j=idx
    %         pfj=bsxfun(@times,refprojs_hat(:,Cjk(k,j),j),shift_phases);
    %         c=real(projs_hat(:,Ckj(k,j),projidx)'*pfj);
    %         bestcorr=max(c);
    %         cvals(j)=bestcorr;
    %     end
    %     score(k)=mean(cvals(idx));
    % end
    % [bestRscore_ref,bestRidx_ref]=max(score);
    % t_ref=toc;
    %
    % log_message('Correlation difference between reference and GPU = %e',bestRscore_ref-bestRscore);
    % log_message('Index difference between reference and GPU = %e',bestRidx_ref-bestRidx);
    % log_message('Correlation difference between reference and CPU = %e',bestRscore_ref-bestRscore_cpu);
    % log_message('Index difference between reference and GPU = %e',bestRidx_ref-bestRidx_cpu);
    % log_message('Timing Ref/CPU/GPU = %5.2f/%5.2f/%5.2f',t_ref,t_cpu,t_gpu);
    %
    % % end of reference code

    if trueRs~=-1
        log_message('\t Frobenius error norm of estimated rotation of projection %d is %5.2e.',...
            projidx,norm(candidate_rots(:,:,bestRidx)-trueRs(:,:,projidx),'fro')/3);
    end
    
    % Now that we have the best shift, the the corresponding 2D shift of the
    % projection.
    shift_equations=zeros(Nrefs,3);
    dtheta=2*pi/L;
    
    idx=find(Mkj(bestRidx,:)==1);
    for j=idx
        pfj=bsxfun(@times,refprojs_hat(:,Cjk(bestRidx,j),j),shift_phases);
        c=real(currproj(:,Ckj(bestRidx,j))'*pfj);
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
if verbose==1
    fprintf('\n');
end
log_message('Total time for orienting %d projections is %5.2f seconds.',...
    n_projs,t_total);
