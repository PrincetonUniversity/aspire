function [sPCA_data, sPCA_coeff, basis, recon_spca ] =  data_sPCA(images, noise_v_r, adaptive_support)
% Tejal April 2016

if nargin < 3 || isempty(adaptive_support)
    adaptive_support = false;
end

n = size(images, 3);
if adaptive_support
    energy_thresh=0.99;
    [ c, R ] = choose_support_v6( cfft2(images), energy_thresh); %Estimate band limit and compact support size
    c=c*(0.5/floor(size(images,1)/2)); % Rescaling between 0 and 0.5
else
    c = 0.5;
    R = floor(size(images, 1)/2);
end
n_r = ceil(4*c*R);
tic_basis=tic;
[ basis, sample_points ] = precomp_fb( n_r, R, c );
timing.basis=toc(tic_basis);
num_pool=10;

log_message('Start computing sPCA coefficients')
[ ~, coeff, mean_coeff, sPCA_coeff, U, D ] = jobscript_FFBsPCA(images, R, noise_v_r, basis, sample_points, num_pool);
log_message('Finished computing sPCA coefficients')

%%The following part is to select top 400 components
Freqs = cell(size(D));
RadFreqs = cell(size(D));
for i = 1:size(D);
    if ~isempty(D{i})
        Freqs{i} = (i-1)*ones(length(D{i}), 1);
        RadFreqs{i} = [1:length(D{i})]';
    end;
end;

Freqs = cell2mat(Freqs);
RadFreqs = cell2mat(RadFreqs);
D = cell2mat(D);
k = min(length(D), 400); %keep the top 400 components
[ D, sorted_id ] = sort(D, 'descend');
D = D(1:k);
Freqs = Freqs(sorted_id(1:k));
RadFreqs = RadFreqs(sorted_id(1:k));
sCoeff = zeros(length(D), n);
for i = 1:length(D)
    sCoeff(i, :) = sPCA_coeff{Freqs(i)+1}(RadFreqs(i), :);
end;

%cumulative energy is 95%, commented out because for the clean data 95% threshold still keeps more than 1000 components. 
%cumD = zeros(size(D, 1), 1);
%for i = 1:length(D)
%    cumD(i) = sum(D(1:i));
%end;
%cumD = cumD/cumD(end);
%threshold = 0.95
%id_cum = find(cumD>threshold, 1)-1;
%length_D = id_cum
%D = D(1:id_cum);
%Freqs = Freqs(1:id_cum);
%RadFreqs = RadFreqs(1:id_cum);

sPCA_data.U = U;
sPCA_data.Freqs = Freqs;
sPCA_data.RadFreqs = RadFreqs;
sPCA_data.Coeff = sCoeff;
sPCA_data.Mean = mean_coeff;
%sPCA_data.FBcoeff = coeff;
sPCA_data.c=c;
sPCA_data.R=R;

L0=size(images,1);
n_max=size(images,3); % Number of images to denoise
%Computes eigen images, need output from IFT_FB.m.
[ fn ] = IFT_FB(R, c);
log_message('Start reconstructing images after sPCA')
[~, recon_spca] = denoise_images_analytical(U, fn, mean_coeff, sPCA_coeff, L0, R, n_max);
log_message('Finished reconstructing images after sPCA')
