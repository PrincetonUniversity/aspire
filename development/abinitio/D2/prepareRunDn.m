
%% initialize some parameters for simulation
max_shift_ratio=0.15;
max_shift=0;%round(129*max_shift_ratio);
shift_step=1;
[projs,Rijs_gt,q,ref_shifts]=genDataForSimulation(vol,...
    100,max_shift,1,1/3,s,10);

%% Generate lookup Data
grid_res=1200;
eq_min_dist=7;
inplane_res=5;
lookup_data=genLookupGrid_eqClass(grid_res,eq_min_dist,inplane_res,s);
[scls_lookup_data]=genSelfCls(lookup_data,2);
[oct1_ij_map,oct2_ij_map]=genSclsScoresIdxMap_eqClass(scls_lookup_data);
scls_lookup_data.oct1_ij_map=oct1_ij_map;
scls_lookup_data.oct2_ij_map=oct2_ij_map;
clear oct1_ij_map oct2_ij_map

%% Run Dn algorithm
doFilter=1;
maxNumWorkers=12;
manyProjs=100;
ntheta=360;
saveDir='/a/home/cc/math/eitanros/cryo/Pipeline/Results/sim';
sampleName='beta_gal_sim_45';
pixA=1.896;
cutoff=0.143;
saveIntermediate=0;

params=struct('max_shift_ratio',max_shift_ratio,'max_shift',max_shift,...
    'shift_step',shift_step,'doFilter',doFilter,'Rijs_gt',Rijs_gt,...,
    'maxNumWorkers',maxNumWorkers,'manyProjs',manyProjs,'ntheta',ntheta,...
    's',s,'q',q,'saveDir',saveDir,'sampleName',sampleName,'vol',vol,...
    'ref_shifts',ref_shifts,'pixA',pixA,'cutoff',cutoff,....
    'saveIntermediate',saveIntermediate,'scl_scores',[],'J_list_in',[]);
stages.st1=0;
stages.st2=0;
stages.st3=1;
stages.st4=1;
stages.st5=0;
stages.Rijs_est=[];
stages.Rijs_synced=[];%results.Rijs_synced;
stages.Rijs_rows=[];%results.Rijs_rows;
stages.colors=[];%results.colors;
stages.rots_est=[];
[results]=runDn(projs,lookup_data,scls_lookup_data,params,stages);
