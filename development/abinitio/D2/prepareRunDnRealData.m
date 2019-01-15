
%% Initialize images and ML data
mapname='/a/home/cc/math/eitanros/cryo/data/averages_nn400_EM_group1.mrcs';
projstack=ReadMRC(mapname);
sidx=1:2000;%9002:2:10000;
projs=projstack(:,:,sidx);
nr=size(projs,1);
max_shift_ratio=0.15;
max_shift=round(nr*max_shift_ratio);
shift_step=1;
vol=ReadMRC(mapname); %Input vol from mrc for comparison
s=rng();

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

%% Initialize parameters for algorithm and run
doFilter=1;
maxNumWorkers=12;
manyProjs=100;
ntheta=360;
saveDir='/a/home/cc/math/eitanros/cryo/Pipeline/Results'; %Directory to save results
sampleName='beta_gal_real_2000_group1'; %Name of molecule for results file
saveIntermediate=0; %Save intermediate stages results
pixA=1.896;
cutoff=0.5;

%parameters for run
params=struct('max_shift_ratio',max_shift_ratio,'max_shift',max_shift,...
    'shift_step',shift_step,'doFilter',doFilter,'Rijs_gt',[],...
    'maxNumWorkers',maxNumWorkers,'manyProjs',manyProjs,'ntheta',ntheta,...
    's',s,'q',[],'saveDir',saveDir,'sampleName',sampleName,'vol',vol,...
    'ref_shifts',[],'pixA',pixA,'cutoff',cutoff,'saveIntermediate',saveIntermediate,...
    'scl_scores',[],'J_list_in',[]);

%Which stages to run
stages.st1=1;
stages.st2=1;
stages.st3=1;
stages.st4=1;
stages.st5=1;

%If some stages are to be skipped, intermediate data should be provided
%Intermediate data is independent, for instance, to begin run from stage 4
%no need to provide data which is needed for stages 2,3.
stages.Rijs_est=[];%results.Rijs_est; %needed for stage 2
stages.Rijs_synced=[];%results.Rijs_synced; %needed for stage 3
stages.Rijs_rows=[];%results.Rijs_rows; %needed for stage 4
stages.colors=[];% results.colors; %also needed for stage 4
stages.rots_est=[]; % needed for stage 5

%Run algorithm
[results]=runDn(projs,lookup_data,scls_lookup_data,params,stages);
