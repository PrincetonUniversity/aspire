function cryo_abinitio_D2(instack,outvol,matfname, gpuIdx,...
    grid_res,eq_min_dist,inplane_res,n_theta,max_shift,...
    shift_step,doFilter)

% Set default parameters

if ~exist('gpuIdx','var') || isempty(gpuIdx)
    nGPU=gpuDeviceCount;
    log_message('No GPU specified. Using all %d GPUs',nGPU);
    gpuIdx=1:nGPU;
end

if ~exist('grid_res','var')
    grid_res=[];
end

if ~exist('eq_min_dist','var')
    eq_min_dist=[];
end

if ~exist('inplane_res','var')
    inplane_res=[];
end

if ~exist('n_theta','var')
    n_theta=[];
end

if ~exist('max_shift','var')
    max_shift=[];
end

if ~exist('shift_step','var')
    shift_step=[];
end

if ~exist('doFilter','var')
    doFilter=[];
end


% Load projections
log_message('Loading mrc stack %s. Plese be patient...', instack);
projs = ReadMRC(instack);
nImages = size(projs,3);
log_message('Done loading mrc image stack file');
log_message('Projections loaded. Using %d projections of size %d x %d',nImages,size(projs,1),size(projs,2));

% Get number of allowed parallel processes
myCluster = parcluster('local');
nCpu = myCluster.NumWorkers;

% Call D2 abinition algorithm
results=runD2(projs,gpuIdx,nCpu,grid_res,eq_min_dist,inplane_res,...
                max_shift,shift_step,n_theta,doFilter);
            
% Save results
log_message('Saving volume %s',outvol);
WriteMRC(results.vol,1,outvol);

log_message('Saving reconstruction parameters %s',matfname);
if ~isempty(matfname)
    save(matfname,'results');
end