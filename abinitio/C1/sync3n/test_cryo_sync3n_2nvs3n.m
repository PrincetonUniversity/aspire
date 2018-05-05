% Compare the accuracy of the 2Nx2N synchronization algorithm with that of
% the 3Nx3N weighted and non-weighted algorithms.
%
% The algorithm generates plots of the mean and median estimation errors.
%
% Yoel Shkolnisky, April 2014.


%% Generate simulated projections
% Generate 200 simulated projections of size 65x65.
% For simplicity, the projections are centered.
initstate;
n=129;
K=1000;
n_theta=360;
n_r=129;
SNRlist=[1/4,1/8,1/16,1/20,1/32,1/48,1/64];
outfile='./test_2nvs3n.mat';

[projs,~,~,rots_ref]=cryo_gen_projections(n,K,100000);
[ref_clstack,~]=clmatrix_cheat(rots_ref,n_theta); % Reference common lines matrix
inv_rots_ref = permute(rots_ref, [2 1 3]);
dirref=R2S2(inv_rots_ref,n_theta);

results=struct;

for snridx=1:numel(SNRlist)
    SNR=SNRlist(snridx);
    results(snridx).SNR=SNR;
    
    noisy_projs=cryo_addnoise(projs,SNR,'gaussian');
    masked_projs=mask_fuzzy(noisy_projs,round(n/2)); % Applly circular mask
    
    
    % Compute polar Fourier transform, using radial resolution n_r and angular
    % resolution n_theta. n_theta is the same as above.
    [npf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
    
    % Find common lines from projections
    max_shift=0;
    shift_step=1;
    clstack = cryo_clmatrix(npf,K,1,max_shift,shift_step);
    prop=comparecl( clstack, ref_clstack, n_theta, 10 );
    fprintf('Percentage of correct common lines: %f%%\n\n',prop*100);
    results(snridx).clrate=prop;
    
    % Estimate orientations using 2Nx2N sychronization.    
    S2N=cryo_syncmatrix_vote(clstack,n_theta);
    rotations2N=cryo_syncrotations(S2N);
    dir2N=R2S2(rotations2N,n_theta);
    [~,~,stats2N]=check_orientations(dirref,dir2N,1,0);
    results(snridx).mean2N=stats2N.mean_err;
    results(snridx).med2N=stats2N.med_err;
    results(snridx).std2N=stats2N.std_err;
    results(snridx).min2N=stats2N.min_err;
    results(snridx).max2N=stats2N.max_err;

    
    % Estimate orientations using non-weighted 3Nx3N synchronization
    rotations3N=cryo_sync3n_estimate_rotations(clstack,n_theta,0);
    dir3N=R2S2(rotations3N,n_theta);
    [~,~,stats3N]=check_orientations(dirref,dir3N,1,0);
    results(snridx).mean3N=stats3N.mean_err;
    results(snridx).med3N=stats3N.med_err;
    results(snridx).std3N=stats3N.std_err;
    results(snridx).min3N=stats3N.min_err;
    results(snridx).max3N=stats3N.max_err;
    
    % Estimate orientations using weighted 3Nx3N synchronization
    rotations3Nw=cryo_sync3n_estimate_rotations(clstack,n_theta,1);
    dir3Nw=R2S2(rotations3Nw,n_theta);
    [~,~,stats3Nw]=check_orientations(dirref,dir3Nw);
    results(snridx).mean3Nw=stats3Nw.mean_err;
    results(snridx).med3Nw=stats3Nw.med_err;
    results(snridx).std3Nw=stats3Nw.std_err;
    results(snridx).min3Nw=stats3Nw.min_err;
    results(snridx).max3Nw=stats3Nw.max_err;

   save(outfile,'results');
   log_message('Saved results for SNR=%5.2e (test %d/%d)',SNR,snridx,numel(SNRlist));
end

%% Plot results
data=load(outfile);
results=data.results;
T=struct2table(results);

% Plot errors
figure;
meandata=table2array(T(:,{'SNR','mean2N','mean3N','mean3Nw'}));
semilogx(meandata(:,1),meandata(:,2),'-rx',meandata(:,1),meandata(:,3),'-bo',meandata(:,1),meandata(:,4),'-gd');
grid on
grid minor
set(gca,'XTick',sort(SNRlist))
xlim([min(SNRlist) max(SNRlist)])
ticLoc = get(gca,'XTick');
ticstr = cellfun(@(x) sprintf('1/%d',round(1/x)),num2cell(ticLoc),'UniformOutput',false);
set(gca,'XTickLabel',ticstr)
set(gca,'Xdir','reverse');
legend('mean 2Nx2N','mean 3Nx3N','mean 3Nx3N weighted');
xlabel('SNR'); ylabel('Error (degrees)');
set(gca,'XTick',sort(SNRlist))

figure;
meddata=table2array(T(:,{'SNR','med2N','med3N','med3Nw'}));
semilogx(meddata(:,1),meddata(:,2),'--rx',meddata(:,1),meddata(:,3),'--bo',meddata(:,1),meddata(:,4),'--gd');
grid on
grid minor
set(gca,'XTick',sort(SNRlist))
xlim([min(SNRlist) max(SNRlist)])
ticLoc = get(gca,'XTick');
ticstr = cellfun(@(x) sprintf('1/%d',round(1/x)),num2cell(ticLoc),'UniformOutput',false);
set(gca,'XTickLabel',ticstr)
set(gca,'Xdir','reverse');
legend('med 2Nx2N','med 3Nx3N','med 3Nx3N weighted','Location','northwest');
xlabel('SNR'); ylabel('Error (degrees)');
