
%% initialize some parameters and data
mapname='/a/home/cc/math/eitanros/cryo/data/averages_nn400_EM_group1.mrcs';
projstack=ReadMRC(mapname);
sidx=1:2000;%9002:2:10000;
projs=projstack(:,:,sidx);
nr=size(projs,1);
masked_projs = projs; %mask_fuzzy(projs,0.5*(nr-1));
max_shift_ratio=0.15;
max_shift=round(nr*max_shift_ratio);
shift_step=1;

% Generate lookup Data
grid_res=1200;
eq_min_dist=15;
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
saveDir='/a/home/cc/math/eitanros/cryo/Pipeline/Results';
sampleName='beta_gal_real_2000_group1';
pixA=1.896;
cutoff=0.5;

params=struct('max_shift_ratio',max_shift_ratio,'max_shift',max_shift,...
    'shift_step',shift_step,'doFilter',doFilter,'Rijs_gt',[],...
    'maxNumWorkers',maxNumWorkers,'manyProjs',manyProjs,'ntheta',ntheta,...
    's',s,'q',[],'saveDir',saveDir,'sampleName',sampleName,'vol',vol,...
    'ref_shifts',[],'pixA',pixA,'cutoff',cutoff,'saveIntermediate',1,...
    'scl_scores',results.scl_scores,'J_list_in',[]);
stages.st1=0;
stages.st2=0;
stages.st3=0;
stages.st4=1;
stages.st5=1;
stages.Rijs_est=[];%results.Rijs_est;
stages.Rijs_synced=[];%results.Rijs_synced;
stages.Rijs_rows=results.Rijs_rows;
stages.colors=results.colors;
stages.rots_est=[];
[results]=runDn(masked_projs,lookup_data,scls_lookup_data,params,stages);
