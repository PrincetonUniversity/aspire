
function [results]=runDn(projs,lookup_data,scls_lookup_data,params,stages)
%% Preapre save dir
if ~isempty(params.saveDir)
    saveDir=params.saveDir;
end
if ~isempty(params.sampleName)
    sampleName=params.sampleName;
end
if ~exist('saveDir','var')
    error('Please input directory to store results');
end

%% Prepare ML data
cls_lookup=[lookup_data.cls_lookup;lookup_data.cls2_lookup];
ntheta=params.ntheta;
nrot=size(projs,3);

%% Prepare shift params
if isempty(params.s)
    rng(params.s);
else
    rng();
end
nr=size(projs,1);
max_shift=params.max_shift;
if isempty(params.max_shift)
    max_shift_ratio=params.max_shift_ratio;
    max_shift=ceil(max_shift_ratio*nr);
end
max_shift_1D = ceil(2*sqrt(2)*max_shift);
shift_step=params.shift_step;

%% Prepare images
masked_projs = mask_fuzzy(projs,0.5*(nr-1));
[pf,~]=cryo_pft(masked_projs,nr,ntheta);

%% Compute self common lines data
doFilter=params.doFilter;
scl_scores_noisy=params.scl_scores;
if isempty(scl_scores_noisy)
    [scl_scores_noisy]=scoreCandidates(scls_lookup_data,pf,max_shift_1D,shift_step,doFilter);
end
scls_lookup_data.scls_scores=scl_scores_noisy;
results.scl_scores=scl_scores_noisy;
if params.saveIntermediate==1
    if exist('sampleName','var')
        save([saveDir,'/',sampleName,'_st0_',cell2mat(datetime_suffix())],'results','-v7.3');
    else
        save([saveDir,'/res_',mat2cell(datatime_suffix)],'results','-v7.3');
    end
end

%% Find first stage to run
run_stages=[stages.st1,stages.st2,stages.st3,stages.st4,stages.st5];
first_stage=find(run_stages,1);

%% Compute common lines data
if stages.st1
    currPool=gcp('nocreate');
    if ~isempty(currPool)
        if currPool.NumWorkers~=4
            delete(gcp('nocreate'));
            gpuDevice([]);
            parpool('local',4);
        end
    end

    [corrs_out]=clmatrix_ML_D2_scls_par(pf,cls_lookup,...
        doFilter,max_shift_1D,shift_step,scls_lookup_data);
    whichScore=2;
    corrs_data.corrs_idx=corrs_out.corrs_idx(:,whichScore);
    corrs_data.ML_corrs=corrs_out.corrs(:,whichScore);
    est_idx=corrs_data.corrs_idx;
    Rijs_est=getRijsFromLinIdx(lookup_data,est_idx);

    results.ML_corrs=corrs_data.ML_corrs;
    results.Rijs_est=Rijs_est;
    results.est_idx=est_idx; %Can be used only with this lookup data

    if ~isempty(params.Rijs_gt)
        Rijs_gt=params.Rijs_gt;
        [sum_err,~,~]=calcApproxErr(Rijs_gt,Rijs_est);
        prc=prctile(sum_err,[55:5:90,91:100]);
        prc2=prctile(sum_err,5:5:50);
        results.prc=prc;
        results.prc2=prc2;
    end
    
    if params.saveIntermediate==1
        if exist('sampleName','var')
            save([saveDir,'/',sampleName,'_st1_',cell2mat(datetime_suffix())],'results','-v7.3');
        else
            save([saveDir,'/res_',mat2cell(datatime_suffix)],'results','-v7.3');
        end
    end
    
elseif first_stage<=2 && stages.st2
    if ~isempty(stages.Rijs_est)
        Rijs_est=stages.Rijs_est;
    else
        error('Cannot continue, missing Rijs_est');
    end
end
%% Run J-synchronization
if stages.st2
    currPool=gcp('nocreate');
    if ~isempty(currPool)
        if currPool.NumWorkers<params.maxNumWorkers && nrot>params.manyProjs
            delete(gcp('nocreate'));
            parpool('local',params.maxNumWorkers);
        end
    end
    J_list_in=params.J_list_in;
    if ~isempty(J_list_in)
        [Rijs_synced,~,evals_jsync]=outer_sync_dev(permute(Rijs_est,[1,2,4,3]),nrot,8,J_list_in);
    else
        [Rijs_synced,~,evals_jsync]=outer_sync_dev(permute(Rijs_est,[1,2,4,3]),nrot,7);
    end
    results.evals_jsync=evals_jsync;
    results.Rijs_synced=Rijs_synced;
    
    if params.saveIntermediate==1
        if exist('sampleName','var')
            save([saveDir,'/',sampleName,'_st2_',cell2mat(datetime_suffix())],'results','-v7.3');
        else
            save([saveDir,'/res_',mat2cell(datatime_suffix)],'results','-v7.3');
        end
    end
elseif first_stage<=3 && stages.st3
    if ~isempty(stages.Rijs_synced)
        Rijs_synced=stages.Rijs_synced;
    else
        error('Cannot continue, missing Rijs_synced');
    end
end

%% Compute coloring
if stages.st3
    currPool=gcp('nocreate');
    if ~isempty(currPool)
        if currPool.NumWorkers<params.maxNumWorkers && nrot>params.manyProjs
            delete(gcp('nocreate'));
            parpool('local',params.maxNumWorkers);
        end
    end 
    %[colors,Rijs_rows,~,D_colors,~,unmix_colors_all]=sync_colors_unmix2(Rijs_synced,nrot,3);
    [colors,Rijs_rows,~,D_colors]=sync_colors_unmix2(Rijs_synced,nrot,3);
    results.D_colors=D_colors;
    results.Rijs_rows=Rijs_rows;
    %results.unmix_colors_all=unmix_colors_all;
    results.colors=colors;
    
    if params.saveIntermediate==1
        if exist('sampleName','var')
            save([saveDir,'/',sampleName,'_st3_',cell2mat(datetime_suffix())],'results','-v7.3');
        else
            save([saveDir,'/res_',mat2cell(datatime_suffix)],'results','-v7.3');
        end
    end
    
elseif first_stage<=4 && stages.st4
    if ~isempty(stages.colors)
        colors=stages.colors;
    else
        error('Cannot continue, missing colors');
    end
    if ~isempty(stages.Rijs_rows)
        Rijs_rows=stages.Rijs_rows;
    else
        error('Cannot continue, missing Rijs_rows');
    end
end

%%  Synchronize signs
if stages.st4
    currPool=gcp('nocreate');
    if ~isempty(currPool)
        if currPool.NumWorkers<params.maxNumWorkers && nrot>params.manyProjs
            delete(gcp('nocreate'));
            parpool('local',params.maxNumWorkers);
        end
    end 
    [rots_est,svals1,svals2,svals3]= restoreRotations_eff_dev2(Rijs_rows,colors(:,1)',nrot);
    results.svals1=svals1;
    results.svals2=svals2;
    results.svals3=svals3;
    results.rots_est=rots_est;
    
    if params.saveIntermediate==1
        if exist('sampleName','var')
            save([saveDir,'/',sampleName,'_st4_',cell2mat(datetime_suffix())],'results','-v7.3');
        else
            save([saveDir,'/res_',mat2cell(datatime_suffix)],'results','-v7.3');
        end
    end
    
else
    if ~isempty(stages.rots_est)
        rots_est=stages.est_rots;
    else
        error('Cannot continue, missing Rijs_rows');
    end
end

%% Analyse results with debug data
q=params.q;
if ~isempty(q)
    debug_params=struct('real_data',0,'analyzeResults',1,'refq',q,'FOUR',4,'alpha',pi/2,'n_theta',ntheta);
    [rot_alligned,err_in_degrees,mse] = analyze_results(rots_est,debug_params);
    results.rot_alligned=rot_alligned;
    results.err_in_degrees=err_in_degrees;
    results.mse=mse;
end

%% Reconstruct volume
if stages.st5
    currPool=gcp('nocreate');
    if ~isempty(currPool) && currPool.NumWorkers<params.maxNumWorkers && nrot>params.manyProjs
        delete(gcp('nocreate'));
        parpool('local',params.maxNumWorkers);
    end
    %[volRec,dxD2] = reconstructDn_ML(noisy_projs,rots,nr,360,max_shift,1);
    [volRec,dxDn,post_corrs] = reconstructDn_dev(masked_projs,rots_est,nr,360,max_shift,1);
    %[volRec,dxD2] = reconstructD2(masked_projs,rots,nr,360,max_shift,1);
    vol=params.vol;
    pixA=params.pixA;
    cutoff=params.cutoff;
    if ~isempty(vol)
        [~,~,volaligned]=cryo_align_densities(vol,volRec,pixA,2);
        results.volaligned=volaligned;
        plotFSC(vol,volaligned,cutoff,pixA);
    end
    results.dxDn=dxDn;
    results.corrs=post_corrs;
end

%% Save results
if exist('saveDir','var') 
    if exist('sampleName','var')
        save([saveDir,'/',sampleName,'_',cell2mat(datetime_suffix())],'results','-v7.3');
    else
        save([saveDir,'/res_',mat2cell(datetime_suffix)],'results','-v7.3');
    end
end


