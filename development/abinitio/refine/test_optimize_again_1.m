% Take projections, estimated rotations, and estimated shifts.

Nprojs=100;
initstate;
q=qrand(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,q);
projs=permute(projs,[2,1,3]);
[projshifted,true_shifts]=cryo_addshifts(projs,[],5,2);
true_shifts=true_shifts.';
snr=100000;
projshifted=cryo_addnoise(projshifted,snr,'gaussian');

% Convert quaternions to rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=(q_to_rot(q(:,k))).';
end

n_r=ceil(size(projs,1)/2);
L=360;
projs_hat=cryo_pft(projshifted,n_r,L,'single');

for k=1:Nprojs
    pf=projs_hat(:,:,k);    
    pf=cryo_raynormalize(pf);
    projs_hat(:,:,k)=pf;
end


dtheta=2*pi/L;
rmax=size(projs_hat,1);
rk2=(0:rmax-1).';
W=zeros(Nprojs*(Nprojs-1)/2,1);
idx=1;

tic;

for k1=1:Nprojs-1
    R1=trueRs(:,:,k1);
    R2=trueRs(:,:,k1+1:Nprojs);
    [Mkj,Ckj,Cjk]=commonline_R_vec(R2,R1,L,0.999);
    
    dx1=true_shifts(:,k1);
      
    for k2=1:Nprojs-k1        
        if Mkj(k2)>0
            
            alpha1=(Cjk(k2)-1)*dtheta;
            delta1=[sin(alpha1) cos(alpha1)]*dx1;
            shift_phases1=exp(-2*pi*sqrt(-1).*rk2*delta1.'./(2*rmax+1));

            U=projs_hat(:,Cjk(k2),k1).*shift_phases1;
            V=projs_hat(:,Ckj(k2),k2+k1);
            V=conj(U).*V;

            dx2=true_shifts(:,k2+k1);
            alpha2=(Ckj(k2)-1)*dtheta;
            delta2=[sin(alpha2) cos(alpha2)]*dx2;
            shift_phases2=exp(-2*pi*sqrt(-1).*rk2*delta2.'./(2*rmax+1));

            V=V.*shift_phases2;
            W(idx)=real(sum(V));

        end
        idx=idx+1;
    end
end
c=sum(W)./sum(W>0);
% Compute the mean correlation using only the valid correlations, that is,
% those for which Mkj>0.

e=1-c;

toc;
disp(e);
