%
% Shift estimatiom for shifted cryo-EM projections.
%
% The is a development code for properly reconstructing a molecule from
% shifted projections.
%
% The function:
%   1. Generates noiseless shifted projections.
%   2. Finds common lines between the projections (taking shifts into
%      account) 
%   3. Generates map2, which contains for each pair of projections (k1,k2)
%      the index of the equation for the relative shift between projections
%      (k1,k2).
%   4. Plots the first singular values of the shift equations. There
%      should be subspace of dimension three with sigular values that are
%      much smaller than the other singular values.
%   5. Solves the shift equations and compares the reference shifts with
%      the estimated ones. The error should small (of the order of 10^-3,
%      due to search in steps of 1 pixel and possible imperfect detection
%      of the correct shift in cryo_clmatrix).
%   6. Reconstructs volumes with and without shift correction, and plots
%      the FSC curves agains the reference volume. The corrected version
%      should have good correlation with the original volume. 
%
% Yoel Shkolnisky, February 2012.

K=200;
L=65;
n_r=100;
n_theta = 360;
[projs,noisy_projs,shifts,rots_ref]=cryo_gen_projections(L,K,1000,6,3);
%[projs,noisy_projs,shifts,rots_ref]=cryo_gen_projections(L,K,1000,0,1);
[pf,sampling_freqs]=cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections   

open_log(0);
[clmatrix,clcorr,shift_equations,shift_equations_map]=cryo_clmatrix(pf,K,1,16,1);

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
est_shifts=shift_equations(:,1:end-1)\shift_equations(:,end);
est_shifts=transpose(reshape(est_shifts,2,K));

s1=reshape(shifts.',2*K,1);
s2=reshape(est_shifts.',2*K,1);
V=V(:,1:end-3); % Null space of shift_equations.

% Compute the difference between the true and estimated shifts in the
% subspace that is orthogonal to the null space of shift_equations.
disp(norm(V.'*(s1-s2))/norm(V.'*s1)); 

inv_rots_ref = permute(rots_ref, [2 1 3]);
dirref=R2S2(inv_rots_ref,n_theta);
% Reconstrut with no shifts
X=cryo_vol2rays(pf);
T1p=2/128;
vol=fastreconstruct3d2_interp(X,dirref,T1p,1,64);

ii=norm(imag(vol(:)))/norm(vol(:));
if ii>1.0e-13
    warning('GCAR:imaginaryComponents','vol hads imaginary components: %7.5e',ii);
end
volr=real(vol);
writeMRC(volr,1,'reshift_test1_a.mrc');

% Reconstruct with shifts
spf=cryo_reshift(pf,sampling_freqs,est_shifts);
X=cryo_vol2rays(spf);
T1p=2/128;
vol=fastreconstruct3d2_interp(X,dirref,T1p,1,64);

ii=norm(imag(vol(:)))/norm(vol(:));
if ii>1.0e-13
    warning('GCAR:imaginaryComponents','vol hads imaginary components: %7.5e',ii);
end
volr=real(vol);
writeMRC(volr,1,'reshift_test1_b.mrc');

% Only reshift_test1_b.mrc looks acceptable.

volref=ReadMRC('cleanrib.mrc');
vol1=ReadMRC('reshift_test1_a.mrc');
vol2=ReadMRC('reshift_test1_b.mrc');


% Registration of vol1 fails.
% dx1=register_translations_3d(vol1,volref,[],1);
% vol1=reshift_vol(vol1,-dx1);
% figure;
% view3d(volref,2.0e-4,'r')
% view3d(vol1,1.0e-4,'g')

% If registering (volref,vol2) then reshift vol2 by dx. If registering
% (vol2,volref) then reshift vol2 by -dx.
dx2=register_translations_3d(vol2,volref,[],1);
vol2=reshift_vol(vol2,-dx2);

figure;
view3d(volref,2.0e-4,'r')
view3d(vol2,2.0e-4,'g')

figure;
plot(FSCorr(vol1,volref),'b');
hold on;
plot(FSCorr(vol2,volref),'r');
hold off;
legend('Not corrected','Corrected');

return;
