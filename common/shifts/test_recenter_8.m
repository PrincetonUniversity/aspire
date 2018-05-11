% Test center of mass correction for noisy projections.
% The script tests how center of mass correction affects the detection rate
% of common lines for noisy images.
% 
% The center of mass is estimated from the noisy projections with and
% without radial masking of the projections.
%
% The projections are NOT centered.
% The script generates a table of the estimation error vs SNR. 
%
% Yoel Shkolsniky, November 2014.

clear;
initstate;

K=100;     % Number of projections.
n=65;     % Size of each projection is nxn.
SNRs=[1000 1 1/2 1/4 1/8 1/16 1/32]; % SNRs to test
n_r=ceil(n/2);  % Radius of each projection
n_theta=360;    % Angular resolution of each projection (for computing its Fourier transform).

% Check if using the following masking in mask_projections
mask_radius=floor(n*0.45);
mask_risetime=floor(n*0.1);
center=(n+1)/2;
mask = fuzzymask([n n],2,mask_radius,mask_risetime,[center center]);
% norm(imag(cfft2(mask))) should be tiny

VERBOSE=0;

fprintf('SNR \t \t NS- \t \t NS+ \t \t CM- \t \t CM+ \t \t MS- \t \t MS+\n');

for j=1:numel(SNRs)
    snr=SNRs(j);
    [projs,noisy_projs,~,rots]=cryo_gen_projections(n,K,snr,3,1);  % Generate projections.
    clstack_ref=clmatrix_cheat(rots,n_theta); % True common lines matrix.

    % Noisy projections with no preprocessing
    pf=cryo_pft(noisy_projs,n_r,n_theta);
    open_log(0);
    
    clstack1=cryo_clmatrix(pf,K,VERBOSE,0,1); % Noisy projections ignore shifts
    prob1=comparecl(clstack1,clstack_ref,n_theta,5);
    clstack2=cryo_clmatrix(pf,K,VERBOSE,10,1); % Noisy projections consider shifts
    prob2=comparecl(clstack2,clstack_ref,n_theta,5);


    % Apply center of mass correction then masking
    projs_aligned=zeros(size(projs));
    cmvec=zeros(K,2);
    for k=1:K
        [projs_aligned(:,:,k),cm]=recenter(noisy_projs(:,:,k));
        projs_aligned(:,:,k)=projs_aligned(:,:,k).*mask;
        cmvec(k,:)=cm;
    end

    pf=cryo_pft(projs_aligned,n_r,n_theta);
    open_log(0);
    
    clstack3=cryo_clmatrix(pf,K,VERBOSE,0,1); % CM corrected ignore shifts
    prob3=comparecl(clstack3,clstack_ref,n_theta,5);
    clstack4=cryo_clmatrix(pf,K,VERBOSE,10,1);% CM corrected consider shifts
    prob4=comparecl(clstack4,clstack_ref,n_theta,5);

    % Apply masking then center of mass correction then masking
    projs_masked=zeros(size(projs));    
    for k=1:K
        masked_proj=noisy_projs(:,:,k).*mask;
        [projs_masked(:,:,k),cm]=recenter(masked_proj);
        projs_masked(:,:,k)=projs_masked(:,:,k).*mask; % Mask again
    end

    pf=cryo_pft(projs_masked,n_r,n_theta);
    open_log(0);

    clstack5=cryo_clmatrix(pf,K,VERBOSE,0,1); % masked+CM corrected ignore shifts
    prob5=comparecl(clstack5,clstack_ref,n_theta,5);
    clstack6=cryo_clmatrix(pf,K,VERBOSE,10,1);% masked+CM corrected consider shifts
    prob6=comparecl(clstack6,clstack_ref,n_theta,5);

%         % Apply filtering then masking then center of mass correction then
%         % masking again
%     projs_filt=zeros(size(projs));    
%     for k=1:K
%         filt_proj=GaussFilt2(noisy_projs(:,:,k),0.2).*mask;
%         [projs_filt(:,:,k),cm]=recenter(filt_proj);
%         projs_filt(:,:,k)=projs_filt(:,:,k).*mask; % Mask again
%     end
% 
%     pf=cryo_pft(projs_filt,n_r,n_theta);
%     open_log(0);
% 
%     clstack7=cryo_clmatrix(pf,K,VERBOSE,10,1); % masked+CM corrected ignore shifts
%     prob7=comparecl(clstack7,clstack_ref,n_theta,5);
    clstack7=clstack6;
    for k=1:K
        % zero all common lines whose center of mass is significantly off
        % center.
        if norm(cmvec(k,:))>10
            clstack_ref(k,:)=0;
            clstack_ref(:,k)=0;
        end
    end
    prob7=comparecl(clstack7,clstack_ref,n_theta,5);
    
    fprintf('%6.4e \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f\n',...
        snr,prob1,prob2,prob3,prob4,prob5,prob6,prob7);

        
end
