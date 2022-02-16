% results_cryo_align_vols_1

clear;
initstate;
% First entry is symmetry group, second is EMD code, third is reported
% model resolution in Angstrom.
test_densities = {...
{'C1',	 2660, 	3.2},...
{'C2',	 0667,	6.2},...
{'C3',	 0731,	2.85},...
{'C4',	 0882,	3.3},...
{'C5',	21376,	2.6},...
{'C7',	11516,	2.38},...
{'C8',	21143,	3.63},...
{'C11',  6458,	4.7},...
{'D2',	30913,	1.93},...
{'D3',	20016,	2.77},...
{'D4',	22462,	2.06},...
{'D7',	 9233,	2.1},...
{'D11',	21140,	3.68},...
{'T',	 4179,	4.1},...
{'O',	22658,	1.36},...
{'I',   24494,	2.27}};
results_1 = cell(numel(test_densities),8);
results_2 = cell(numel(test_densities),8);
for testidx = 1:numel(test_densities)   
    %% Generate two density maps. 
    test_data = test_densities{testidx};
    symmetry = test_data{1};
    emdid = test_data{2};
    res = test_data{3};
    log_message('********************************');
    log_message('Test %d/%d %s  (EMD%04d)',testidx,numel(test_densities),symmetry,emdid);  
    try
        mapfile = cryo_fetch_emdID(emdid);
        vol = ReadMRC(mapfile);
        vol = GaussFilt(vol,0.1);
        sz_vol = 64;
        vol = cryo_downsample(vol,sz_vol,0);
    catch err
        delete(mapfile);
    end
    delete(mapfile);
    [R,~,~] = svd(rand(3)); % Generate random rotation
    volRotated1 = fastrotate3d(vol,R); % Rotate the reference volume by the random rotation
    
    volRotated2 = flip(volRotated1,3);
    volRotated2 = reshift_vol(volRotated2,[-5 0 0]);
    
    volRotated1 = reshift_vol(volRotated1,[-5 0 0]);   
    
    %% Align volumes:
    G = genSymGroup(symmetry);
    n = size(G,3);
    opt.N_projs = 30;     
    opt.G = G;
    opt.true_R = R;
    opt.dofscplot = 0;
    opt.downsample = 64;
    opt.sym = symmetry;

    tic;
    [Rest,estdx1,reflect1,vol2aligned,bestcorr1,T_optimize] = cryo_align_vols(vol,volRotated1,1,opt);
    T = toc;
    c2 = corr(vol(:),vol2aligned(:));
    %% Accurate error calculation:
    % The difference between the estimated and reference rotation should be an
    % element from the symmetry group:
    n_g = size(G,3);
    g_est_t = R.'*Rest.';
    dist = zeros(n_g,1);
    for g_idx = 1:n_g
        dist(g_idx,1) = norm(g_est_t-G(:,:,g_idx),'fro');
    end
    [~,min_idx] = min(dist);
    g_est = G(:,:,min_idx);
    vec_ref = rotationMatrixToVector(R.');
    angle_ref = norm(vec_ref);
    axis_ref = vec_ref/angle_ref;    
    vec_est = rotationMatrixToVector(g_est*Rest);
    angle_est = norm(vec_est);
    axis_est = vec_est/angle_est;  
    Err_ang1 = acosd(dot(axis_est,axis_ref)); %deg
    Err_ang2 = abs(rad2deg(angle_ref)-rad2deg(angle_est)); %deg

    results_1{testidx,1} = emdid;
    results_1{testidx,2} = symmetry; 
    results_1{testidx,3} = reflect1;
    results_1{testidx,4} = c2; %Correlation between the aligned volumes
    results_1{testidx,5} = T; %Time of alignment
    results_1{testidx,6} = T_optimize; %Time of optimization
    results_1{testidx,7} = Err_ang1; %Angle between axes
    results_1{testidx,8} = Err_ang2; %Angle difference
    
    
    %% Align volumes with reflection:
    tic;
    [Rest,estdx,reflect,vol2aligned,bestcorr,T_optimize] = cryo_align_vols(vol,volRotated2,1,opt);
    T = toc;
    c2 = corr(vol(:),vol2aligned(:));
    %% Accurate error calculation:
    % The difference between the estimated and reference rotation should be an
    % element from the symmetry group:
    n_g = size(G,3);
    g_est_t = R.'*Rest.';
    dist = zeros(n_g,1);
    for g_idx = 1:n_g
        dist(g_idx,1) = norm(g_est_t-G(:,:,g_idx),'fro');
    end
    [~,min_idx] = min(dist);
    g_est = G(:,:,min_idx);
    vec_ref = rotationMatrixToVector(R.');
    angle_ref = norm(vec_ref);
    axis_ref = vec_ref/angle_ref;    
    vec_est = rotationMatrixToVector(g_est*Rest);
    angle_est = norm(vec_est);
    axis_est = vec_est/angle_est;  
    Err_ang1 = acosd(dot(axis_est,axis_ref)); %deg
    Err_ang2 = abs(rad2deg(angle_ref)-rad2deg(angle_est)); %deg

    results_2{testidx,1} = emdid;
    results_2{testidx,2} = symmetry; 
    results_2{testidx,3} = reflect;
    results_2{testidx,4} = c2; %Correlation between the aligned volumes
    results_2{testidx,5} = T; %Time of alignment
    results_2{testidx,6} = T_optimize; %Time of optimization
    results_2{testidx,7} = Err_ang1; %Angle between axes
    results_2{testidx,8} = Err_ang2; %Angle difference  
    
end    
    
    