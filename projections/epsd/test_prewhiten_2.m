% Test the performance of the PSD projection method in whitening colored noise.

N = 32;
K = 1000;

initstate;

% Generate a stack of noise images
noise1 = noise_exp2d(N, K, 0);
noise2 = 1/sqrt(2)*noise_rexpr(N, K, 0);

% Geneate random coefficients.
a1 = randn(K, 1);
a2 = randn(K, 1);

% Combine the noise sources using the random coefficients to get the actual
% noise images.
noise = bsxfun(@times, noise1, permute(a1, [2 3 1])) + ...
    bsxfun(@times, noise2, permute(a2, [2 3 1]));

% Calculate periodograms for all the images.
per = cryo_periodogram(noise, 1:N^2);

% Get mean and covariance.
per_mean = cryo_periodogram_mean(per);
per_cov = cryo_periodogram_covariance(per);

% Extract top eigenvectors and project onto the subspace spanned.
[V, ~] = mdim_eigs(per_cov, 2, 'la');
per_proj = cryo_project_psd(per, V_oracle);
per_proj = max(per_proj, 0);

% Whiten using the mean PSD estimate and the individual PSD estimates obtained
% from projection.
whitened_mean = cryo_prewhiten(noise, per_mean, 5e-2);
whitened_proj = cryo_prewhiten(noise, per_proj, 5e-2);

% Calculate periodogram of the whitened noise images to see how well we did.
per_whitened_mean = cryo_periodogram(whitened_mean, 1:N^2);
per_whitened_proj = cryo_periodogram(whitened_proj, 1:N^2);

% Compare the periodogram obtained from the whitened images with the goal:
% a PSD of all ones. Note that this error will never be zero due to the
% measurement noise inherent in the periodogram. The lowest possible value is
% 1 for the complex coefficients and 2 for the real ones (i.e. frequencies 0 
% and N/2).
mse_mean = mean(abs(per_whitened_mean-ones([N*ones(1, 2) K])).^2, 3);
mse_proj = mean(abs(per_whitened_proj-ones([N*ones(1, 2) K])).^2, 3);

% Display the resulting MSE profiles.
figure();
plot([0:N/2-1]/N, mse_mean(N/2+1:end,N/2+1), ...
    [0:N/2-1]/N, mse_proj(N/2+1:end,N/2+1));
legend('MSE for mean periodogram', 'MSE for projected periodogram');
