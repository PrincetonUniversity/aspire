% Print a table of detection rates using cryo_clmatrix_cpu
%
% Yoel Shkolnisky, June 2016.

fprintf('Commonline detection rates.\n');
fprintf('===========================\n');
n = 89;
K=100;
n_r=100;
n_theta=360;
for k=1:6
    SNR=1/2^k;
    fprintf('\nSNR=1/%d, ',2^k)
    max_shift2d=0;
    step_size=1;
    initstate;
    [p, np, shifts, rots] = ...
        cryo_gen_projections(n,K,SNR,max_shift2d,step_size);
    [ref_clmatrix,clcorr]=clmatrix_cheat(rots,n_theta);
    

    mask_radius=round(size(np,1)*0.45);
    [np,sigma]=mask_fuzzy(np,mask_radius);
    [npf,freqs]=cryo_pft(np,n_r,n_theta);
    max_shift1d = ceil(2*sqrt(2)*max_shift2d);
    shift_step = 1;
    [ clstack,corrstack, shift_equations,shift_equations_map]...
                        = cryo_clmatrix_cpu( npf,K,1,max_shift1d,shift_step,10 );
    prop=comparecl( clstack, ref_clmatrix, n_theta, 10 );
    fprintf('correctness of common lines: %f\n\n',prop);
    fprintf('-----------------------------------------------------------------\n');
end

