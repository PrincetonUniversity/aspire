% SHRINK_COVAR Shrink eigenvalues of a covariance matrix
%
% Usage
%    B = shrink_covar(A, noise_var, gamma, shrinker);
%
% Input
%    A: The sample covariance matrix to be shrunk.
%    noise_var: The variance of the noise in the original samples.
%    gamma: The gamma parameter determining p/n, where p is the dimension
%       of the space and n is the number of samples.
%    shrinker: Specifies the shrinker type (default 'frobenius_norm'):
%          - 'operator_norm': The operator norm shrinkage proposed by Donoho
%             et al [1].
%          - 'frobenius_norm': The Frobenius norm shrinkage proposed by Donoho
%             et al. [1].
%          - 'soft_threshold': A soft thresholding shrinkage, where lambda is
%             mapped to lambda-(1+sqrt(gamma))^2.
%          - 'none': Do not perform any shrinkage.
% Output
%    B: The covariance matrix `A` with its eigenvalues shrunk according to
%       `noise_var`, `gamma` and `shrinker`.
%
% References:
%    [1] D. L. Donoho, M. Gavish, and I. M. Johnstone, Optimal shrinkage of
%        eigenvalues in the spiked covariance model, arXiv preprint
%        arXiv:1311.0851, (2013).

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function B = shrink_covar(A, noise_var, gamma, shrinker)
    if nargin < 4 || isempty(shrinker)
        shrinker = 'frobenius_norm';
    end

    shrinkers = {'operator_norm', 'frobenius_norm', 'soft_threshold', ...
        'none'};

    if ~ismember(shrinker, shrinkers)
        error('Invalid `shrinker`.');
    end

    A = A/noise_var;

    [V, D] = eig(make_symmat(A));

    lambda_max = (1+sqrt(gamma))^2;

    D(D<lambda_max) = 0;

    if strcmp(shrinker, 'operator_norm')
        lambda = D(D>lambda_max);
        lambda = 1/2*(lambda-gamma+1+sqrt((lambda-gamma+1).^2-4*lambda))-1;
        D(D>lambda_max) = lambda;
    elseif strcmp(shrinker, 'frobenius_norm')
        lambda = D(D>lambda_max);
        lambda = 1/2*(lambda-gamma+1+sqrt((lambda-gamma+1).^2-4*lambda))-1;
        c = (1-gamma./lambda.^2)./(1+gamma./lambda);
        lambda = lambda.*c;
        D(D>lambda_max) = lambda;
    elseif strcmp(shrinker, 'soft_threshold')
        lambda = D(D>lambda_max);
        D(D>lambda_max) = lambda-lambda_max;
    end

    B = V*D*V';

    B = B*noise_var;
end
