
K=250;
L=65;
n_r=100;
n_theta = 360;
[projs,noisy_projs,shifts,rots_ref]=cryo_gen_projections(L,K,3,6,3);
%
% Corrupt two projections
projs(10:50,10:50,1:5)=0;
noisy_projs(10:50,10:50,1:5)=0;
masked_projs=mask_fuzzy(noisy_projs,35);
[cpf,sampling_freqs]=cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections   
[npf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections   

open_log(0);
% Find common lines using noisy projections, but reoncstruct using clean.
[clmatrix,clcorr,shift_equations,shift_equations_map]=cryo_clmatrix(npf,K,1,16,1);

[PHIc,eigs_diary,removed_projections,clean_clmatrix]=...
    cryo_clmat2orientations(clmatrix,zeros(size(clmatrix)),n_theta,...
    size(clmatrix,1),10,'/tmp/yoel');

save tmp;


[FE,ES,SF]=filter_shift_equations(shift_equations,removed_projections,clean_clmatrix,shifts);

map2=zeros(K);
idx=1;
for k1=1:K  
    for k2=k1+1:K
        map2(k1,k2)=idx;
        idx=idx+1;
    end
end

if norm(map2-shift_equations_map)>1.0e-15
    error('wrong map');
end

[U,S,V]=svd(full(shift_equations(:,1:end-1)));
s=diag(S); % The first three singular values should be much smaller than the others.
figure;
bar(s(end:-1:end-19))


% Center projections
 clidx=find(clean_clmatrix~=0);  % common-lines present in the final common-lines matrix.
 [r,c]=ind2sub(size(clean_clmatrix),clidx);
 clidx=clidx(c>r);
 [I,J]=ind2sub(size(clean_clmatrix),clidx);
 
 retained_projections=setdiff(1:K,removed_projections);
 retained_projections=retained_projections(:);
 retained_cl=sub2ind([K K],retained_projections(I),retained_projections(J));
 retained_shifts=sort(map2(retained_cl));
 proj_idx=([2*(retained_projections-1)+1 2*retained_projections]).';
 proj_idx=[proj_idx(:) ; 2*K+1];
 
 % Now retained_shifts contained the indices of the relevant common-line
 % equations. proj_idx are the indices of the shift variables to retain.
 shift_equations=shift_equations(retained_shifts,proj_idx);

[U,S,V]=svd(full(shift_equations(:,1:end-1)));
s=diag(S); % The first three singular values should be much smaller than the others.
figure;
bar(s(end:-1:end-19))
  
 est_shifts=shift_equations(:,1:end-1)\shift_equations(:,end);
 est_shifts=transpose(reshape(est_shifts,2,length(retained_projections)));

 % Filter the reference shifts.
 idx=proj_idx(1:end-1);
 shifts_filtered=shifts.';
 shifts_filtered=shifts_filtered(:);
 shifts_filtered=shifts_filtered(idx);
 shifts_filtered=transpose(reshape(shifts_filtered,2,length(retained_projections)));
% est_shifts2=shift_equations(:,1:end-1)\shift_equations(:,end);
% est_shifts2=transpose(reshape(est_shifts,2,K));

s1=reshape(shifts_filtered.',2*numel(retained_projections),1);
s2=reshape(est_shifts.',2*numel(retained_projections),1);
V=V(:,1:end-3); % Null space of shift_equations.

% Compute the difference between the true and estimated shifts in the
% subspace that is orthogonal to the null space of shift_equations.
disp(norm(V.'*(s1-s2))/norm(V.'*s1)); 

disp(norm(FE(:)-shift_equations(:)));
disp(norm(SF(:)-shifts_filtered(:)));
disp(norm(ES(:)-est_shifts(:)));

cpf(:,:,removed_projections)=[];
npf(:,:,removed_projections)=[];
rots_ref(:,:,removed_projections)=[];
inv_rots_ref = permute(rots_ref, [2 1 3]);
dirref=R2S2(inv_rots_ref,n_theta);
PHIc=register_orientations(PHIc,dirref);
% Reconstrut with no shifts
X=cryo_vol2rays(cpf);
T1p=2/128;
vol=fastreconstruct3d2_interp(X,PHIc,T1p,1,64);

ii=norm(imag(vol(:)))/norm(vol(:));
if ii>1.0e-13
    warning('GCAR:imaginaryComponents','vol hads imaginary components: %7.5e',ii);
end
volr=real(vol);
writeMRC(volr,1,'reshift_test2_a.mrc');

% Reconstruct with shifts
spf=cryo_reshift(cpf,sampling_freqs,est_shifts);
X=cryo_vol2rays(spf);
T1p=2/128;
vol=fastreconstruct3d2_interp(X,PHIc,T1p,1,64);

ii=norm(imag(vol(:)))/norm(vol(:));
if ii>1.0e-13
    warning('GCAR:imaginaryComponents','vol hads imaginary components: %7.5e',ii);
end
volr=real(vol);
writeMRC(volr,1,'reshift_test2_b.mrc');

% Only reshift_test1_b.mrc looks acceptable.

volref=ReadMRC('cleanrib.mrc');
vol1=ReadMRC('reshift_test2_a.mrc');
vol2=ReadMRC('reshift_test2_b.mrc');

dx1=register_translations(vol1,volref);
vol1=reshift_vol(vol1,dx1);
dx2=register_translations(vol2,volref);
vol2=reshift_vol(vol2,dx2);

figure;
plot(FSCorr(vol1,volref),'b');
hold on;
plot(FSCorr(vol2,volref),'r');
hold off;
legend('Not corrected','Corrected');

return;
