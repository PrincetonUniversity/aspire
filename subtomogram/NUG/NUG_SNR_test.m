cd /u/liuyuan/Documents/MATLAB/Subtomo_average/NUG_cryoEM_forYuan
useNUG_yuan

%% projections to Fourier space
nSNR = 8; rep = 5;
MSE = zeros(3,nSNR, rep);
T = zeros(3,nSNR, rep);

%% generate subtomogram tilt series
n = 64;
ns = 5;
K_sub = 8;
K = K_sub*ns;

% The quaternions are in the form q = cos(t/2) + (x,y,z)sin(t/2).
b_sub = pi/24;
initstate
subq = qrand(K_sub);
halfns = (ns-1)/2;
tiltq = zeros(4,K);   
grid = repmat((-halfns:halfns)*b_sub,1,K_sub);     % tilt angles in radians, their relative
                                        % angle is 7.5 deg = 1/24*pi
tiltq(1,:) = cos(grid/2);   %%%%%% use grid/2 to fit the quaternion definition            
tiltq(2,:) = sin(grid/2);         % assume the rotation axis of tilt series 
                                % is the x-axis, the projection is in the z-axis
ref_q = zeros(halfns,K); % refq = b*a = q*r which is rotation b follows rotation a in absolute frame ### check
for i = 1:K
    a = subq(:,ceil(i/ns));
    b = tiltq(:,i);
    ref_q(1,i) = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4);
    ref_q(2,i) = a(1)*b(2) + a(2)*b(1) - a(3)*b(4) + a(4)*b(3);
    ref_q(3,i) = a(1)*b(3) + a(2)*b(4) + a(3)*b(1) - a(4)*b(2);
    ref_q(4,i) = a(1)*b(4) - a(2)*b(3) + a(3)*b(2) + a(4)*b(1);
end
% 
% % Generate clean projections
% [projs,~,~,refq]=cryo_gen_projections(n,K,10^4,[],[],[],ref_q,[]);
% masked_projs=mask_fuzzy(projs,33); % Applly circular mask
% viewstack(masked_projs,ns,K_sub);
rots_ref = q_to_rot(ref_q);
rots_sub = q_to_rot(subq);
%%

for iSNR = 1:nSNR
    for j = 1:rep
        SNR = 10^4/10^(iSNR);
        %SNR = 100;
        [~,noisy_projections,~,refq]=cryo_gen_projections(n,K,SNR,[],[],[],rots_ref,[]);

        % cryo NUG
        tic
        [MSE(1,iSNR, j),~] = cryoEM_blackBox_y(noisy_projections, rots_ref);
        T(1,iSNR, j) = toc;

        tic
        [MSE(2,iSNR, j),~] = subtomo_blackBox(noisy_projections, rots_sub, K_sub,ns,b_sub);
        T(2,iSNR, j) = toc;

        n_theta = 72;
        n_r = n/2;
        tic
        % Calculate polar Fourier transforms of images.
        pf = cryo_pft(noisy_projections, n_r, n_theta);

        % Estimate common-lines matrix.
        clstack_est = cryo_clmatrix_cpu(pf);

        % Construct syncronization matrix from common lines.
        S_est = cryo_syncmatrix_vote(clstack_est, n_theta);

        % Estimate rotations using synchronization.
        [~,~,MSE(3,iSNR, j),~] = cryo_syncrotations(S_est, rots_ref);
        T(3,iSNR, j) = toc;
        
%         % using subtomogram for sychronization
%         [clstack, ~, ~, ~] = commonlines_gaussian_vsub(npf,max_shift,shift_step);


        save(sprintf('nugSNR2019_%dim_rep5.mat',K));
    end
end

%%
% K = 40;
% load(sprintf('nugSNR2019_%dim_rep5.mat',K));
% SNR = 10^4./10.^(1:nSNR);
% figure(1)
% subplot(2,1,1)
% semilogx(SNR./n^2,mean(MSE(1,:,:),2),'x-'); hold on;
% semilogx(SNR./n^2,mean(MSE(2,:,:),2),'x-'); hold on;
% semilogx(SNR./n^2,mean(MSE(3,:,:),2),'x-'); hold on;
% legend('Single Projection with NUG','Subtomograms with NUG','Synchronization');
% xlabel('Normalized SNR');
% ylabel('MSE');
% title({'MSE for estimated orientation with N=40 images of size 64'});
% subplot(2,1,2)
% semilogx(SNR./n^2,mean(T(1,:,:),2),'x-'); hold on;
% semilogx(SNR./n^2,mean(T(2,:,:),2),'x-'); hold on;
% semilogx(SNR./n^2,mean(T(3,:,:),2),'x-'); hold on;
% legend('Single Projection','Subtomograms','Synchronization');
% xlabel('Normalized SNR');
% ylabel('Time (s)');
% title({'Time to estimate orientation with 40 images of size 64, on 12 cores'});
% %set(gca,'FontSize',16)
% % h=gcf;
% % set(h,'PaperOrientation','landscape');
% print('-fillpage','nug_snr','-dpdf')


