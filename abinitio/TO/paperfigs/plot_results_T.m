%% Results for T symmetry
% Yoel Shkolnisky, November 2022

[~,s]=ReadMRC('./results_T/volref.mrc');
pixA = s.pixA;
cutoff = 0.5;

% 25 images, variable SNR
vols={'volref.mrc','vol_25_1_aligned.mrc',...
    'volref.mrc','vol_25_2_aligned.mrc',...
    'volref.mrc','vol_25_3_aligned.mrc',...
    'volref.mrc','vol_25_4_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'1000','1','1/2','1/4'});
figname=fullfile(results_dir,'fsc_T_N25.eps');
print('-depsc',figname);

% 50 images, variable SNR
vols={'volref.mrc','vol_50_1_aligned.mrc',...
    'volref.mrc','vol_50_2_aligned.mrc',...
    'volref.mrc','vol_50_3_aligned.mrc',...
    'volref.mrc','vol_50_4_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'1000','1','1/2','1/4'});
figname=fullfile(results_dir,'fsc_T_N50.eps');
print('-depsc',figname);


% 100 images, variable SNR
vols={'volref.mrc','vol_100_1_aligned.mrc',...
    'volref.mrc','vol_100_2_aligned.mrc',...
    'volref.mrc','vol_100_3_aligned.mrc',...
    'volref.mrc','vol_100_4_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'1000','1','1/2','1/4'});
figname=fullfile(results_dir,'fsc_T_N100.eps');
print('-depsc',figname);

% 200 images, variable SNR
vols={'volref.mrc','vol_200_1_aligned.mrc',...
    'volref.mrc','vol_200_2_aligned.mrc',...
    'volref.mrc','vol_200_3_aligned.mrc',...
    'volref.mrc','vol_200_4_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'1000','1','1/2','1/4'});
figname=fullfile(results_dir,'fsc_T_N200.eps');
print('-depsc',figname);


% snr_list = [1000, 1, 1/2, 1/4, 1/8];
% SNR = 1000; variable number of images
vols={'volref.mrc','vol_25_1_aligned.mrc',...
    'volref.mrc','vol_50_1_aligned.mrc',...
    'volref.mrc','vol_100_1_aligned.mrc',...
    'volref.mrc','vol_200_1_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'25','50','100','200'});
figname=fullfile(results_dir,'fsc_T_SNR1.eps');
print('-depsc',figname);

% SNR = 1; variable number of images
vols={'volref.mrc','vol_25_2_aligned.mrc',...
      'volref.mrc','vol_50_2_aligned.mrc',...
    'volref.mrc','vol_100_2_aligned.mrc',...
    'volref.mrc','vol_200_2_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'25','50','100','200'});
figname=fullfile(results_dir,'fsc_T_SNR2.eps');
print('-depsc',figname);

% SNR = 1/2; variable number of images
vols={'volref.mrc','vol_25_3_aligned.mrc',...
    'volref.mrc','vol_50_3_aligned.mrc',...
    'volref.mrc','vol_100_3_aligned.mrc',...
    'volref.mrc','vol_200_3_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'25','50','100','200'});
figname=fullfile(results_dir,'fsc_T_SNR3.eps');
print('-depsc',figname);


% SNR = 1/4; variable number of images
vols={'volref.mrc','vol_25_4_aligned.mrc',...
    'volref.mrc','vol_50_4_aligned.mrc',...
    'volref.mrc','vol_100_4_aligned.mrc',...
    'volref.mrc','vol_200_4_aligned.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end
plotFSCmany(vols,cutoff,pixA,{'25','50','100','200'});
figname=fullfile(results_dir,'fsc_T_SNR4.eps');
print('-depsc',figname);

% Print some noisy images
clf;
for snr=1:4
    mrcName = fullfile(results_dir,sprintf('noisy_projs_25_%d.mrcs',snr));
    projs=ReadMRC(mrcName);
    for n=1:8
        imagesc(projs(:,:,n));
        colormap(gray)
        axis image
        axis off
        projs_name = fullfile(results_dir,sprintf('noisy_T_n%d_snr%d.eps',n,snr));
        print('-depsc',projs_name);
    end
end
        
% Plot timing
timing=zeros(4,1);
N=[25,50,100,200];
for k=1:numel(N)
    total_time = 0.0;
    for snr=1:4
        fid=fopen(fullfile(results_dir,sprintf('timing_%d_%d.txt',N(k),snr)),'r');
        t = fscanf(fid,'%g');
        total_time = total_time + t;
        fclose(fid);
    end
    timing(k)=total_time/4;
    disp(timing(k))
end
figure;
plot(N,timing,'x-r','LineWidth',2);
xticks(N)
xlim([min(N)-5,max(N)+5])
grid on
%coeffs = polyfit(N(:),timing(:),2);
%hold on;
%plot(N,polyval(coeffs,N))
%hold off;

% Plot FSC of relion's initial model
mapfile = cryo_fetch_emdID(10835); % Download 3D map file with O symmetry
vol_ref=ReadMRC(mapfile);           % Load the map file
pixA=0.639 ;
results_dir = './results_T';
WriteMRC(vol_ref,pixA,fullfile(results_dir,'volref_10389.mrc'));
vol_relion = ReadMRC('/data/yoelsh/datasets/10389/relion/InitialModel/job058/run_it300_class001.mrc');
[~,~,vol_relion_aligned]=cryo_align_densities(vol_ref,vol_relion,pixA,1,0.5);
WriteMRC(vol_relion_aligned,pixA,fullfile(results_dir,'vol_relion_aligned_10389.mrc'));
vol_aspire = ReadMRC('/data/yoelsh/datasets/10389/abinitio/T_symmetry_data_set_10389_with_2084_candidates.mrc');
[~,~,vol_aspire_aligned]=cryo_align_densities(vol_ref,vol_aspire,pixA,1,0.5);
WriteMRC(vol_aspire_aligned,pixA,fullfile(results_dir,'vol_aspire_aligned_10389.mrc'));

vols={'volref_10389.mrc','vol_aspire_aligned_10389.mrc',...
    'volref_10389.mrc','vol_relion_aligned_10389.mrc'};

results_dir = './results_T';
for k=1:length(vols)
    vols{k}=fullfile(results_dir,vols{k});
end

plotFSCmany(vols,0.5,pixA,{'aspire','relion'});
print('-depsc',fullfile(results_dir,'fsc_10389.eps'));

%plotFSC(vol_ref,vol_relion_aligned,0.5,pixA);
%print('-depsc',fullfile(results_dir,'fsc_abinitio_relion_10389.eps'));

%plotFSC(vol_ref,vol_aspire_aligned,0.5,pixA);
%print('-depsc',fullfile(results_dir,'fsc_abinitio_aspire_10389.eps'));


% figure;
% imagesc(squeeze(sum(GaussFilt(vol_ref,0.05),3)))
% colormap(gray)
% axis image
% axis off
% print('-depsc','ref_view.eps');
% 
% figure;
% imagesc(squeeze(sum(GaussFilt(vol_relion_aligned,1.0),3)));
% colormap(gray)
% axis image
% axis off
% print('-depsc','relion_initial_model_view.eps');
% 
% 
% figure;
% imagesc(squeeze(sum(GaussFilt(vol_aspire_aligned,1.0),3)));
% colormap(gray)
% axis image
% axis off
% 
