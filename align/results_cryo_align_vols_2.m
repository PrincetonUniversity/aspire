% results_cryo_align_vols_2 
%
% Align two density maps where one is rotate/translated relative to
% the other, in the presence of noise.

clear;
initstate;

%% Generate two density maps.

load cleanrib
mapfile = cryo_fetch_emdID(2660); %C1 
%mapfile = cryo_fetch_emdID(0744); %C1 reported resolution=2.82
vol = ReadMRC(mapfile);vol = GaussFilt(vol,0.8);

[R,~,~] = svd(rand(3)); % Generate random rotation
volRotated = fastrotate3d(vol,R); % Rotate the reference volume by the random rotation
volRotated = reshift_vol(volRotated,(rand(3,1)-1/2)*8); % Shift +/- 4 pixels.

vol = vol/norm(vol(:));
volRotated = volRotated/norm(volRotated(:));

SNRs = [100 1 1/2 1/8 1/32 1/128 1/256 1/512].';
corrsNoisy = zeros(numel(SNRs),1);
corrsClean = zeros(numel(SNRs),1);

results = cell(numel(SNRs),7);

for snridx = 1:numel(SNRs)
    
    snr = SNRs(snridx);
    sigma = sqrt(1/snr);
    n1 = sigma*randn(size(vol))./sqrt(numel(vol(:)));
    n2 = sigma*randn(size(volRotated))./sqrt(numel(volRotated(:)));
    volNoisy = vol+n1;
    volRotatedNoisy = volRotated+n2;
           
    %% Visulaize noisy level
    figure
    imagesc(volNoisy(:,:,floor((size(volNoisy,3)+1)/2))); 
    colormap(gray);
    axis image;
    axis off;
    
    %% Align 
    initstate;
    sz_vol = 64;
    vol_copy = cryo_downsample(vol,sz_vol,0);
    volRotated_copy = cryo_downsample(volRotated,sz_vol,0);
    volRotatedNoisy = cryo_downsample(volRotatedNoisy,sz_vol,0);
    
    verbose=1;
    opts.true_R=R;
    opts.sym = 'C1';
    opts.downsample = 64;
    opts.N_projs = 30; 
    
    %[Rest,estdx,volAlignedNoisy]=cryo_align_densities(volNoisy,volRotatedNoisy,0,verbose,0.5,R,0);
    tic;
    [Rest,estdx,reflect,volAlignedNoisy,bestcorr,T_optimize] = cryo_align_vols(vol_copy,volRotatedNoisy,verbose,opts);
    T = toc;
    % Correlation of noisy volumes after alignment
    c1 = corr(vol_copy(:),volAlignedNoisy(:));
    fprintf('Correlation between original and aligned noisy volumes %7.4f\n',c1);
    corrsNoisy(snridx) = c1;
    
    % Correlation of clean volumes according to estimated alignment
    % parameters
    
    volAligned = fastrotate3d(volRotated_copy,Rest);
    volAligned = reshift_vol(volAligned,estdx);
    
    c2=corr(vol_copy(:),volAligned(:));
    fprintf('Correlation between original and aligned clean volumes %7.4f\n',c2);
    corrsClean(snridx) = c2;
    
    % The original and aligned volume should have high agreement according to
    % the FSC. FSC may degrade due to sub-pixel misalignment errors.
    plotFSC(vol_copy,volAlignedNoisy,0.5,1);
    plotFSC(vol_copy,volAligned,0.5,1);
    
    %% Accurate error calculation:
    % The difference between the estimated and reference rotation should be an
    % element from the symmetry group:
    G = genSymGroup('C1');
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

    results{snridx,1} = 1/snr; %1/SNR level
    results{snridx,2} = c1; %Correlation between the noisy aligned volumes
    results{snridx,3} = c2; %Correlation between the clean aligned volumes
    results{snridx,4} = T; %Time of alignment
    results{snridx,5} = T_optimize; %Time of optimization
    results{snridx,6} = Err_ang1; %Angle between axes
    results{snridx,7} = Err_ang2; %Angle difference
    
    
end
%fprintf('1/SNR  corrNoisy  corrClean\n');
%disp([1./SNRs corrsNoisy corrsClean]);


