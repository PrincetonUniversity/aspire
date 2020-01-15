function [results] = runD2(projs,gpuIdx,nCpu,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,ntheta,doFilter,s,Rijs_gt,q,debugParam)
%  runD2 computes a reconstruction of a D2 volume from a given set of
%  projection images.

%  Input parameters: 
%
%  projs       a stack of 2D projection images (3 x 3 x #(projections)).
%  gpuIdx      indices of gpu devices to use for common line detection.
%  nCpu        size of matlab pool for parallel CPU cores. Default value is 1  
%  grid_res     number of sampling points on sphere for projetion directions. 
%               These are generated using the Saaf - Kuijlaars algoithm.
%               Default value is 1200. 
%  inplane_res   The sampling resolution of in-plane rotations for each
%                projetion direction. Default value is 5. 
%  eq_min_dist   Width of strip around equator projection directions from
%                which we DO NOT sample directions. Default value is 7. 
%  ntheta       number of fourier rays on the polar FT of each 2D
%               image. Default value is 360. 
%  max_shift    max 2D shift to try when estimating common lines. Default
%               value is computed as 15% of the projection image dimensions. 
%  shift_step   resolution of shifts to try when estimating common lines.
%               Default value is 1. 
%  doFilter     Use a gaussian mask to filter 2D fourier transforms of 
%               projection when estimating common lines. Default value is
%               1 (do the filtering). 
%  Rijs_gt      Ground truth relative rotations. Use for debuging in
%               simulations. Optional. 
%  s            A pre generated seed ( s = rng() ) for reproducibility.
%               Optional. 
%  q            In a simulation this is an array of unit quaternions from 
%               which to retrieve the rotations used to simulate the 
%               projections. Optional. 
%  debugParam   parameters used to compare resulting volume with a given 
%               volume (using Fourier shell correlations). Optional. 
            
%% Initialize default values for unspecified arguments. 
if ~exist('grid_res','var')
    grid_res = 1200;
end

if ~exist('eq_min_dist','var')
    eq_min_dist = 7;
end

if ~exist('inplane_res','var')
    inplane_res = 5;
end

if ~exist('shift_step','var')
    shift_step = 1200;
end

if ~exist('ntheta','var')
    ntheta = 360;
end

if ~exist('doFilter','var')
    doFilter = 1;
end

if ~exist('nCpu','var')
   nCpu = 1; 
end

if exist('s','var')
    rng(s);
else
    s = rng();
    results.s = s;
end

%% Generate lookup Data
%  genSclsScoresIdxMap_eqClass generates a lookup data for maximum
%  likelihood like scheme to estimate common lines between the images. 
lookup_data=genLookupGrid_eqClass(grid_res,eq_min_dist,inplane_res);
[scls_lookup_data]=genSelfCls(lookup_data,2);
cls_lookup=[lookup_data.cls_lookup;lookup_data.cls2_lookup];
nrot=size(projs,3);

%% Prepare shift params
nr=size(projs,1);
if ~exist('max_shift','var')
    max_shift=ceil(0.15*nr);
end
max_shift_1D = ceil(2*sqrt(2)*max_shift);

%% Prepare images
masked_projs = mask_fuzzy(projs,0.5*(nr-1));
[pf,~]=cryo_pft(masked_projs,nr,ntheta);

%% Compute self common lines data
[scl_scores_noisy]=...
    cryo_clmatrix_scl_ML(pf,max_shift_1D,shift_step,doFilter,scls_lookup_data);
scls_lookup_data.scls_scores=scl_scores_noisy;
results.scl_scores=scl_scores_noisy;

%% Estimate relative rotations using ML procedure on common lines. 
[corrs_out]=clmatrix_ML_D2_scls_par(pf,cls_lookup,...
    doFilter,max_shift_1D,shift_step,scls_lookup_data,gpuIdx);
corrs_data.corrs_idx=corrs_out.corrs_idx;
corrs_data.ML_corrs=corrs_out.corrs;
est_idx=corrs_data.corrs_idx;
Rijs_est=getRijsFromLinIdx(lookup_data,est_idx);

results.ML_corrs=corrs_data.ML_corrs;
results.Rijs_est=Rijs_est;
results.est_idx=est_idx; %Can be used only with this lookup data

if  exist('Rijs_gt','var')
    [sum_err,~,~]=calcApproxErr(Rijs_gt,Rijs_est);
    prc=prctile(sum_err,[55:5:90,91:100]);
    prc2=prctile(sum_err,5:5:50);
    results.prc=prc;
    results.prc2=prc2;
end

%% Create new pool after we are finished with GPU's. 
if nCpu > 1
    currPool = parpool('local',nCpu);
    nWorkers = currPool.NumWorkers;
else
    nWorkers = 1;
end
%% Run J-synchronization (Handedness synchronization).    
[Rijs_synced,~,evals_jsync]=jsync(permute(Rijs_est,[1,2,4,3]),nrot,nWorkers);
results.evals_jsync=evals_jsync;
results.Rijs_synced=Rijs_synced;

%% Synchronize colors
[colors,Rijs_rows,~,D_colors]=syncColors(Rijs_synced,nrot,nWorkers);
results.D_colors=D_colors;
results.Rijs_rows=Rijs_rows;
results.colors=colors;
   
%%  Synchronize signs
[rots_est,svals1,svals2,svals3]= syncSigns(Rijs_rows,colors(:,1)',nrot);
results.svals1=svals1;
results.svals2=svals2;
results.svals3=svals3;
results.rots_est=rots_est;
    
%% Analyse results with debug data
if exist('q','var')
    debug_params=struct('real_data',0,'analyzeResults',1,'refq',q,'FOUR',4,'alpha',pi/2,'n_theta',ntheta);
    [rot_alligned,err_in_degrees,mse] = analyze_results(rots_est,debug_params);
    results.rot_alligned=rot_alligned;
    results.err_in_degrees=err_in_degrees;
    results.mse=mse;
end

%% Reconstruct volume
[volRec,dxDn,post_corrs] = reconstructDn_dev(masked_projs,rots_est,nr,360,max_shift,1);
results.vol=volRec;
results.dxDn=dxDn;
results.corrs=post_corrs;

%% Compare with existing volume. 
if exist('debugParam','var')
    vol=debugParam.vol;
    pixA=debugParam.pixA;
    cutoff=debugParam.cutoff;
    
    [~,~,volaligned]=cryo_align_densities(vol,volRec,pixA,2);
    results.volaligned=volaligned;
    plotFSC(vol,volaligned,cutoff,pixA);
end

