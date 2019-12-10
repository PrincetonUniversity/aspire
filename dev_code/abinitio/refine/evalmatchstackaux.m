function e=evalmatchstackaux(X,projs_hat,L)


Nprojs=size(projs_hat,3);
% Unpack X into rotations and translations
Rs=zeros(3,3,Nprojs); % Each row is (psi,theta,phi) of one of the rotations
eulerangs=reshape(X(1:3*Nprojs),Nprojs,3);
for k=1:Nprojs
    psi=eulerangs(k,1);
    theta=eulerangs(k,2);
    phi=eulerangs(k,3);
    Rs(:,:,k)=Rz(phi)*Ry(theta)*Rx(psi);
end

%norm(Rs(:)-refR(:))

dxs=zeros(2,Nprojs);
for k=1:Nprojs
    dxs(1,k)=X(3*Nprojs+2*(k-1)+1);
    dxs(2,k)=X(3*Nprojs+2*(k-1)+2);
end

%norm(dxs(:)-refdx(:))

dtheta=2*pi/L;
rmax=size(projs_hat,1);
rk2=(0:rmax-1).';
W=zeros(Nprojs*(Nprojs-1)/2,1);
idx=1;

tic;

for k1=1:Nprojs-1
    R1=Rs(:,:,k1);
    R2=Rs(:,:,k1+1:Nprojs);
    [Mkj,Ckj,Cjk]=commonline_R_vec(R2,R1,L,0.999);
    
    dx1=dxs(:,k1);
    alpha1=(Cjk-1)*dtheta;
    delta1=[sin(alpha1) cos(alpha1)]*dx1;
    shift_phases1=exp(-2*pi*sqrt(-1).*rk2*delta1.'./(2*rmax+1));    
    for k2=1:Nprojs-k1        
        if Mkj(k2)>0            
            U=projs_hat(:,Cjk(k2),k1).*shift_phases1(:,k2);
            V=projs_hat(:,Ckj(k2),k2+k1);
            V=conj(U).*V;

            dx2=dxs(:,k2+k1);
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
