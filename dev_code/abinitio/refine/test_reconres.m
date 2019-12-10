% Check the effect of noise on the resolution of the reconstruction.
%
% Generate simulated projection images at various levels of noise,
% reconstruct from the noisy projections and the true orientations (no
% shifts are added) and compare the resolution of the reconstructed volume
% to that of the original (clean) one.
%
% Yoel Shkolnisky, August 2017.

Nprojs=50000;
initstate;
rots = rand_rots(Nprojs);  % Generate Nprojs projections to orient.
voldata=load('cleanrib');
projs=cryo_project(voldata.volref,rots);
projs=permute(projs,[2,1,3]);
% Invert rotations
trueRs=zeros(3,3,Nprojs);
for k=1:Nprojs
    trueRs(:,:,k)=rots(:,:,k).';
end
n=size(projs,1);


snrlist=[10^6 1 1/2 1/4 1/8 1/16 1/32 1/164 1/128 1/256];
FSCs=zeros(floor(n/2),numel(snrlist));
for k=1:numel(snrlist)
    log_message('Processing SNR %d/%d',k,numel(snrlist));
    snr=snrlist(k);
    noisy_projs=cryo_addnoise(projs,snr,'gaussian');
    
    [ v1, ~, ~ ,~, ~, ~] = recon3d_firm( noisy_projs,trueRs, [],1e-8, 100, zeros(n,n,n));
    ii1=norm(imag(v1(:)))/norm(v1(:));
    log_message('Relative norm of imaginary components = %e\n',ii1);
    v1=real(v1);
    FSCs(:,k)=FSCorr(voldata.volref,v1);
    log_message('resolution=%7.4f',fscres(FSCs(:,k),0.143));
    
    X=1:size(FSCs,1);
    plot(repmat(X(:),1,numel(snrlist)),FSCs)
    f1 = @(k) sprintf('   %.3f (%5.3fA)',snrlist(k),fscres(FSCs(:,k),0.143));
    legend(cellfun(f1, num2cell(1:numel(snrlist)), 'UniformOutput', false))
    print(sprintf('res%d.eps',Nprojs),'-depsc');
end

