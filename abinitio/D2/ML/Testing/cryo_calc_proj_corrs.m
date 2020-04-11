
function [C,est_shifts]=cryo_calc_proj_corrs(pf,doFilter,max_shift,shift_step,B)

T=size(pf,2);
if mod(T,2)~=0
    error('n_theta must be even');
end

%band limit
if nargin>4
    pf=pf(1:B+1,:,:);
end

%pf=[flipdim(pf(2:end,T/2+1:end,:),1) ; conj(pf(:,1:T/2,:)) ];
n_theta=size(pf,2);
n_proj=size(pf,3);

rmax=(size(pf,1)-1)/2;
rk=-rmax:rmax; rk=rk(:);
H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2));
H=repmat(H(:),1,n_theta);  % Filter for common-line detection.

% Bandpass filter and normalize each ray of each projection.
% XXX We do not override pf since it is used to debugging plots below. Once
% XXX these debugging plots are removed, replace pf3 by pf. This will save
% XXX a lot of memory.
pf3=pf;
if doFilter
    for k=1:n_proj
        proj=pf(:,:,k);
        proj=proj.*H;
        proj(rmax:rmax+2,:)=0;
        proj=cryo_raynormalize(proj);
        pf3(:,:,k)=proj;
    end
else
    for k=1:n_proj
        proj=pf(:,:,k);
        proj(rmax:rmax+2,:)=0;
        proj=cryo_raynormalize(proj); % QUESTION: Do we need this in the clean case???        
        pf3(:,:,k)=proj;
    end
end


proj1=pf3(:,:,1);
P1=proj1(1:rmax,:);  % Take half ray plus the DC
P1_flipped=conj(P1);

%P1=cat(2,P1,P1_flipped);

n_shifts=ceil(2*max_shift/shift_step+1);
rk=-rmax:rmax; rk=rk(:);
rk2=rk(1:rmax);
% Prepare the shift_phases
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end

P1_stack=zeros(size(P1,1),2*size(P1,2)*n_shifts);
for k=1:n_shifts
    P1_stack(:,(2*k-2)*n_theta+1:(2*k-1)*n_theta)=bsxfun(@times,P1,shift_phases(:,k));
    %P1_flipped_stack(:,(k-1)*n_theta+1:k*n_theta)=bsxfun(@times,P1_flipped,shift_phases(:,k));
    P1_stack(:,(2*k-1)*n_theta+1:(2*k)*n_theta)=bsxfun(@times,P1_flipped,shift_phases(:,k));
end

if norm(proj1(rmax+1,:))>1.0e-13
    error('DC component of projection is not zero');
end

proj2=pf3(:,:,2); % proj1 and proj2 are both normalized to unit norm.
P2=proj2(1:rmax,:);
%P2_flipped=conj(P2);
%P2=cat(2,P2,P2_flipped);
%g_P2=gpuArray(single(P2));

if norm(proj2(rmax+1,:))>1.0e-13
    error('DC component of projection is not zero');
end

C=2*real(P1_stack'*P2);
C=reshape(C,2*n_theta,n_shifts,n_theta);
C=permute(C,[1,3,2]);
[C,est_shifts]=max(C,[],3);
C=squeeze(C);
est_shifts=-max_shift+(squeeze(est_shifts)-1)*shift_step;

%C=[2*real(P1'*P2),2*real(P1_flipped'*P2)];




