function [ class, class_refl, rot, corr, timing ] = Initial_classification(data, r_max, n_nbor, isrann )
%This function does initial classfication on preprocessed data.
%   Input: 
%           data:  data should be phase flipped, prewhitened and roughly centered.
%           r_max: the radius of the region of interest.
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
% Zhizhen Zhao Updated 2013

% Get date dimension:
L=size(data, 1);
P=size(data, 3);
N=floor(L/2);
trunc = 150; % maximum number of FBsPCA leading componnents chosen. 
  
[x, y]=meshgrid(-N:N, -N:N);    %generate square grids
r=sqrt(x.^2+y.^2);      %compute the radius at each grid point
test=reshape(data, L^2, P);
test=test(r>r_max, :);
noise_variance=var(test(:)); %Estimate the noise variance from ear region.

%Use Fourier Bessel steerable PCA to denoise and compress the original data
%set.
tic_FBsPCA=tic;
[ U, D, freqs, rad_freqs, Mean ] = FB_SPCA(data, r_max);
[ UU, Freqs, ~, W ] = FBSPCA_MP_rankEst(P, U, D, freqs, rad_freqs, max(noise_variance, D(trunc)));
toc_FBsPCA=toc(tic_FBsPCA);
%Cutoff is determined by the maximum of the noise_variance and the 200th
%eigenvalues.
fprintf('\nThe number of components is %d \n', size(Freqs, 1));

tic_WF = tic;
% Using Linear Wiener Filter to denoise:
[ Coeff ] = WF_FBSPCA( data, Mean, r_max, [0, 0], UU, Freqs, W, 1);

Coeff(Freqs==0, :)=Coeff(Freqs==0, :)/sqrt(2);
for i=1:P
    Coeff(:, i)=Coeff(:, i)/norm(Coeff(:, i));
end;
Coeff(Freqs==0, :)=Coeff(Freqs==0, :)*sqrt(2);
toc_WF=toc(tic_WF);

[ Coeff_b, Coeff_b_r, toc_bispec ] = Bispec_2Drot( Coeff, Freqs );

if P<=10000
    %%For small dataset, search for nearest neighbors
    tic_nn=tic;
    corr=real(Coeff_b'*[Coeff_b, Coeff_b_r]);
    corr=corr-sparse(1:P, 1:P, ones(P, 1), P, 2*P);
    [~, class]=sort(corr(1:P, :), 2, 'descend');
    class=class(:, 1:n_nbor);
    toc_nn=toc(tic_nn);
    %This part is the brute force nearest neighbor search.
else
    tic_nn=tic;
    if ~isrann
        %Brute Force
        P_max=2000;
        for i=1:ceil(P/P_max)
            corr=real(Coeff_b(:, (i-1)*P_max+1:min(i*P_max, P))'*[Coeff_b, Coeff_b_r]);
            [~, tmp]=sort(corr, 2, 'descend');
            class((i-1)*P_max+1:min(i*P_max, P), :) = tmp(:, 2:n_nbor+1);
            fprintf('\nCorrelation step %d\n', i);
        end;
    else
        %Randomized approximate nearest neighbor search
        [class, ~]=test_points_from_file64([Coeff_b, Coeff_b_r], n_nbor, 10, 0, 0);
    end;
    toc_nn=toc(tic_nn);
end

% Rotational alignment for nearest neighbor pairs
k_max=max(Freqs);
Cell_Coeff=cell(k_max+1, 1);
for i=1:k_max+1
    Cell_Coeff{i}=[Coeff(Freqs==i-1, :), conj(Coeff(Freqs==i-1, :))];
end;
list=[class(:), repmat([1:P]', n_nbor, 1)];

tic_rot=tic;
[corr, rot ]=rot_align(Freqs, Cell_Coeff, list);
toc_rot=toc(tic_rot);

corr=reshape(corr, P, n_nbor);
rot=reshape(rot, P, n_nbor);
class=reshape(class, P, n_nbor);
[corr, id_corr]=sort(corr, 2, 'descend');

for i=1:P
    class(i, :)=class(i, id_corr(i, :));
    rot(i, :)=rot(i, id_corr(i, :));
end;

class_refl=ceil(class/P);
class(class>P)=class(class>P)-P;

timing.FBsPCA=toc_FBsPCA;
timing.WF = toc_WF;
timing.bispec=toc_bispec;
timing.nn=toc_nn;
timing.rot=toc_rot;



