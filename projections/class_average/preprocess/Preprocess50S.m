% Preprocess the 50S dataset

%% Load the data
projs=ReadImagic('~/data/work/molec/datasets/RibosomeData/Rib_flipcor_1.img');
%projs=projections;
sz=size(projs);
n=sz(1);
K=sz(3);

%% Rough  centering of the images

% Normalize - subtract mean of each image and normalize to energy 1.
parfor k=1:K
    p=projs(:,:,k);
    m=mean(p(:));
    p=p-m;
    p=p./norm(p(:));
    projs(:,:,k)=p;
end

transSD=100;
ccprior=-Radius(n).^2/(2*transSD^2);

rmask=fuzzymask(n,2,0.4*n,0.1*n);  % basic mask
for i=1:3  % 3 rounds of centering
    ref=CircSymmetrize(mean(projs,3)).*rmask;
    subplot(1,3,1);
    imagesc(ref);
    axis image;
    title(['Iteration ' num2str(i)]);
    drawnow;
    [projs, shifts, amps]=TransAlignInt(projs,ref,ccprior);
    shiftsr=sqrt(sum(shifts.^2));
%     ok=(shiftsr'<10 & amps>26 & amps<70);  % hard-coded selection
    %     sum(ok)
    subplot(1,3,2);
    plot(shifts(1,:),shifts(2,:),'k.');
    axis([-n/2 n/2 -n/2 n/2]);
    axis equal;
    title(i);
    subplot(1,3,3);
    plot(amps);
    drawnow;
end;

%% Images must be of odd size

if mod(n,2)==0
    projs=projs(1:n-1,1:n-1,:);
    n=n-1;
end

% Normalize again
projs0=zeros(n,n,K);
parfor k=1:K
    p=projs(:,:,k);
    m=mean(p(:));
    p=p-m;
    p=p./norm(p(:));
    projs0(:,:,k)=p;
end


%% Prewhiten the images
mean_data=mean(projs0, 3);
data=projs0;
parfor k=1:K
    data(:, :, k)=data(:, :, k)-mean_data;
end;
[ P ] = Noise_Estimation( data );
[ P2,negmask2 ] = radial_average( P );
[ data ] = cryo_prewhiten(data, P2);

figure;
subplot(1,2,1);
plot(P2(n+1,:));
title('1D profile of the noise spectrum');
subplot(1,2,2);
imagesc(log(-negmask2)); colorbar;
title('Negative (zeroed) entries in the noise spectrum');


% Bring back the mean data
% parfor k=1:K
%     data(:, :, k)=data(:, :, k)+mean_data;
% end;

% Verify that the power spectrum is flat
data0=data;
mean_data0=mean(data0, 3);
parfor k=1:K
    data0(:, :, k)=data0(:, :, k)-mean_data0;
end;

[ P3 ] = Noise_Estimation( data0 );
[ P3,negmask3 ] = radial_average( P3 );

figure;
subplot(1,2,1);
plot(P3(n+1,:));
title('Verification - 1D profile of the noise spectrum');
subplot(1,2,2);
imagesc(log(-negmask3)); colorbar;
title('Verification - Negative (zeroed) entries in the noise spectrum');


% Mask all projections
rmask=fuzzymask(n,2,0.4*n,0.1*n);  % basic mask
projs1=zeros(size(projs));
parfor k=1:K
    projs1(:,:,k)=data(:,:,k).*rmask;
end
    
K=2000;
pf=cryo_pft(projs1(:,:,1:K),50,360);
max_shift=4;
shift_step=1;
tic;
[ clstack,corrstack, shift_equations,shift_equations_map3]=cryo_clmatrix_gpu(pf,K,1,max_shift,shift_step);
toc

return;
