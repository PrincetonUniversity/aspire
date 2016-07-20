function [ class, class_refl, rot, corr,  timing ] = Initial_classification_FD(sPCA_data, n_nbor, isrann )
%Description:
%This function does initial classfication on preprocessed data.
%   Input: 
%           sPCA: steerable PCA data.
%           n_nbor: number of nearest neighbors used for class averaging.
%           isrann: using randomized algorithm for finding k nearest
%           neighbors. isrann=0: brute force nearest neighbor search.
%           isrann=1: randomized nearest neighbor search
%   Output:
%           class: a matrix of size n x n_nbor. This provides a list of
%           nearest neighbors for each image.
%           rot: rotational alginment
%           corr: normalized cross correlation between nearest neighbors
%           timing: timing of different steps in the algorithm
%	    FBsPCA_data: FBsPCA denoising data. It includes r_max--radius of 
%           region of interest, UU--eigenvectors, Freqs--associated angular 
%           frequence, Coeff--image expansion coefficients on UU, Mean--estimated
%           rotationally invariant mean of the data, and W--weight for wiener type filtering.
% Zhizhen Zhao Updated Jan 2015


Coeff = [sPCA_data.Coeff, conj(sPCA_data.Coeff)]; % Tejal April 2016
Freqs = sPCA_data.Freqs;
%rad_Freqs = sPCA_data.rad_Freqs;
n_im = (size(Coeff, 2))/2; % Tejal April 21, 2017
%normalize the coefficients
Coeff(Freqs==0, :)=Coeff(Freqs==0, :)/sqrt(2);
for i=1:2*n_im  %% Tejal April 21, 2017
    Coeff(:, i)=Coeff(:, i)/norm(Coeff(:, i));
end;
Coeff(Freqs==0, :)=Coeff(Freqs==0, :)*sqrt(2);


%Compute bispectrum
%[ Coeff_b, toc_bispec ] = Bispec_2Drot_large( Coeff, Freqs ); %If the number of images and number of Coefficients are large use Bispec_2Drot_large
[ Coeff_b,  toc_bispec ] = Bispec_2Drot_1( Coeff, Freqs );
sprintf('Tejal Size of coeff_b is %d', size(Coeff_b))

if n_im<=10000
    %%For small dataset, search for nearest neighbors
    tic_nn=tic;
    corr=real(Coeff_b(:, 1:n_im)'*Coeff_b);
    %corr=real((Coeff_b(:, 1:n_im))'*[Coeff_b, Coeff_b_r]); % Tejal April 21 2016
    corr=corr-sparse(1:n_im, 1:n_im, ones(n_im, 1), n_im, 2*n_im);
    [~, class]=sort(corr(1:n_im, :), 2, 'descend');
    class=class(:, 1:n_nbor);
    toc_nn=toc(tic_nn);
    %This part is the brute force nearest neighbor search.
else
    tic_nn=tic;
    if ~isrann
        %Brute Force
        P_max=2000;
        for i=1:ceil(n_im/P_max)
            corr = real(Coeff_b(:, (i-1)*P_max+1:min(i*P_max, n_im))'*Coeff_b);
            %corr=real(Coeff_b(:, (i-1)*P_max+1:min(i*P_max, P))'*[Coeff_b, Coeff_b_r]); % Tejal April 21 2016
            [~, tmp]=sort(corr, 2, 'descend');
            class((i-1)*P_max+1:min(i*P_max, n_im), :) = tmp(:, 2:n_nbor+1);
            fprintf('\nCorrelation step %d\n', i);
        end;
    else
        %Randomized approximate nearest neighbor search
        [class, ~]=test_points_from_file64(Coeff_b, n_nbor, 10, 0, 0);
    end;
    toc_nn=toc(tic_nn);
end

% Rotational alignment for nearest neighbor pairs
k_max=max(Freqs);
Cell_Coeff=cell(k_max+1, 1);
for i=1:k_max+1
    Cell_Coeff{i}= Coeff(Freqs==i-1, :);
end;
list=[class(:), repmat([1:n_im]', n_nbor, 1)];

%Initial in-plane rotational alignment within nearest neighbors
tic_rot=tic;
[corr, rot ] = rot_align(max(Freqs), Cell_Coeff, list);
toc_rot=toc(tic_rot);

corr = reshape(corr, n_im, n_nbor);
rot=reshape(rot, n_im, n_nbor);
class=reshape(class, n_im, n_nbor);
[corr, id_corr]=sort(corr, 2, 'descend');

for i=1:n_im
    class(i, :)=class(i, id_corr(i, :));
    rot(i, :)=rot(i, id_corr(i, :));
end;

class_refl=ceil(class/n_im);
class(class>n_im)=class(class>n_im)-n_im;

timing.bispec=toc_bispec;
timing.nn=toc_nn;
timing.rot=toc_rot;

