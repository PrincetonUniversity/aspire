% Basic script to test the performance of the PSD projection method when
% different noise spectra are present.

N = 32;
Ks = [100 1000 10000];

initstate;

% Generate component fixed noise images.
noise1 = noise_exp2d(N, max(Ks), 0);
noise2 = 1/sqrt(2)*noise_rexpr(N, max(Ks), 0);

% Generate random coefficients.
a1 = randn(max(Ks), 1);
a2 = randn(max(Ks), 1);

% Calculate periodograms of the fixed noises.
per1 = cryo_periodogram(noise1, 1:N^2);
per2 = cryo_periodogram(noise2, 1:N^2);

% This is not exact, but in the absence of a proper power spectrum reference,
% it will have to do.
Sref1 = cryo_periodogram_mean(per1);
Sref2 = cryo_periodogram_mean(per2);

% Combine the noise sources using the random coefficients to get the actual
% noise images.
noise = bsxfun(@times, noise1, permute(a1, [2 3 1])) + ...
    bsxfun(@times, noise2, permute(a2, [2 3 1]));

% Calculate the periodograms of the noise images.
per = cryo_periodogram(noise, 1:N^2);

for j = 1:numel(Ks)
    fprintf('Processing K=%d\n', Ks(j));

    % Estimate PSD mean and covariance from periodograms.
    per_mean = cryo_periodogram_mean(per(:,:,1:Ks(j)));
    per_cov = cryo_periodogram_covariance(per(:,:,1:Ks(j)));

    % Extract the top 2 eigenvectors of the covariance and project PSD
    % estimates.
    [V, ~] = mdim_eigs(per_cov, 2, 'la');
    per_proj = cryo_project_psd(per(:,:,1:Ks(j)), V);

    % Compute "oracle" subspace as the span of the fixed power spectra and
    % project.
    V_oracle = cat(3, Sref1, Sref2);
    per_oracle = cryo_project_psd(per(:,:,1:Ks(j)), V_oracle);

    % Compute the true power spectra for each image.
    Sref = bsxfun(@times, Sref1, permute(abs(a1(1:Ks(j))).^2, [2 3 1])) + ...
        bsxfun(@times, Sref2, permute(abs(a2(1:Ks(j))).^2, [2 3 1]));

    % Baseline estimator: just assign the mean PSD to all images.
    per_mean = repmat(per_mean, [1 1 Ks(j)]);

    % Compute the relative 2-norm of the estimates.
    err_mean(j) = norm(per_mean(:)-Sref(:))/norm(Sref(:));
    err_proj(j) = norm(per_proj(:)-Sref(:))/norm(Sref(:));
    err_oracle(j) = norm(per_oracle(:)-Sref(:))/norm(Sref(:));
end
