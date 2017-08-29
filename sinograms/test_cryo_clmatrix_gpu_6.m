% Print a table of detection rates using cryo_clmatrix_gpu
%
% The function demonstrates the effectiveness of correlation map averaging
% for high levels of noise.
%
% Yoel Shkolnisky, June 2016.
clear;
initstate;
fprintf('Commonline detection rates.\n');
fprintf('===========================\n');
K=100;
n_r=100;
n_theta=360;
n=129;
SNRlist=[1/32 1/48 1/64];
for k=1:numel(SNRlist)
    SNR=SNRlist(k);
    fprintf('\nSNR=1/%d, ',1/SNR)
    max_shift_2d=0;
    shift_step_2d=1;
    
     %data=load('cleanrib.mat');
     data.volref=ReadMRC('/home/idog/matlab/molecules/frank70S_g1_nn50_n129.mrc');
     volref=real(data.volref);
     volref=cryo_downsample(volref,[n,n,n]);
     initstate;
     rots_ref = rand_rots(K);
     projs=cryo_project(volref,rots_ref,n);
     projs=permute(projs,[2 1 3]); % Swap dimensions for compitability with old gen_projections.
     [projs,~]=cryo_addshifts(projs,[],max_shift_2d,shift_step_2d);
     noisy_projs=cryo_addnoise(projs,SNR,'gaussian');

    
%     [p, np, shifts, rots_ref] = ...
%         cryo_gen_projections(n,K,SNR,max_shift2d,step_size_2d);
    [ref_clmatrix,clcorr]=clmatrix_cheat(rots_ref,n_theta);
    

    mask_radius=round(size(noisy_projs,1)*0.45);
    [noisy_projs,sigma]=mask_fuzzy(noisy_projs,mask_radius);
    [npf,freqs]=cryo_pft(noisy_projs,n_r,n_theta);
    max_shift1d = ceil(2*sqrt(2)*max_shift_2d);
    shift_step = 1;
    
    % No averaging of the correlation map;
    [ clstack1,corrstack1, shift_equations1,shift_equations_map1]...
                        = cryo_clmatrix( npf,K,1,max_shift1d,shift_step);
    prop1=comparecl( clstack1, ref_clmatrix, n_theta, 10 );
    
    % With averaging of the correlation map.
    filter_radius2=5;
    [ clstack2,corrstack2, shift_equations2,shift_equations_map2]...
                        = cryo_clmatrix( npf,K,1,max_shift1d,shift_step,filter_radius2);
    prop2=comparecl( clstack2, ref_clmatrix, n_theta, 10 );
    
    % With averaging of the correlation map3
    filter_radius3=10;
    [ clstack3,corrstack3, shift_equations3,shift_equations_map3]...
                        = cryo_clmatrix( npf,K,1,max_shift1d,shift_step,filter_radius3);                    
    prop3=comparecl( clstack3, ref_clmatrix, n_theta, 10 );

    fprintf('correctness of common lines: %f (no averaging)  %f (with averaging 5) %f (with averaging 10)\n\n',prop1,prop2,prop3);
    fprintf('-----------------------------------------------------------------\n');
end

