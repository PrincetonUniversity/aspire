
function [P, sigma, Rsquare, Pij, scores_hist, scores_fit, cum_scores] =...
    triangles_scores(Rij, scores_hist, Pmin, Pmax, verbose)
% Triangles Scores
% [P, sigma, Rsquare, Pij, scores_hist, cum_scores, plotted] =
%   triangles_scores(Rij, N, required_P, PLOT, verbose)
% 
% This is a wrapper for the mex functions triangles_scores_mex() &
% pairs_probabilities_mex().
% It mainly runs the functions and analyzes their output to deduce some
% indicators for the reliability of the relative rotations data.
% More information can be found in the mex functions documentations.
% 
% Input:
% Rij - relative rotations matrices (3 X 3 X N-choose-2)
% N - number of images
% required_P - minimum required value for estimated rate of indicative common lines.
% If value is lower, then the function returns in the middle (to save time) and the
% caller should improve the data and call again. Set to 0 to allow any P.
% PLOT - whether to plot the histogram & its fit.
% 
% Output:
% P - estimated rate of indicative common lines
% sigma - typical error in indicative common lines estimations [degrees]
% Rsquare - R^2 value of the fit of the histogram which derives P & sigma
% Pij - estimated probabilities of the common lines to be indicative
% scores_hist - histogram of the common lines triangles scores
% cum_scores - pairs cumulated triangles scores
% plotted - whether the histogram was plotted or not (depends on both PLOT
% and the flow of the function)
% 
% Written by Ido Greenberg, 2015

% Input parsing
if ~exist('Pmin','var'); Pmin=0; end
if ~exist('Pmax','var'); Pmax=1; end
if ~exist('verbose','var'); verbose=true; end
if ~exist('scores_hist','var'); scores_hist=[]; end

% Histogram Analysis Constants
hist_intervals = 100; % number of intervals in the scores histogram
a = 2.2; % empirical constant factor in the theoretical formula of the scores distribution
% NOTE: original analysis found 'a' to equal 2.2. other variants of the code used 2.0 and 2.5.
% The value is supposed to be global, though it might slightly depend on N
% due to biases of the voting algorithm. 2.5 seemed optimal at some
% experiment of N~100.
peak2sigma = 2.43e-2; % empirical relation between the location of the peak of the histigram, and the mean error in the common lines estimations

% Initialization
N=(1+sqrt(1+8*size(Rij,3)))/2; % Extract N from the number of relative rotations.
    % The number of relative rotations for a given N is (N choose 2).
    % Solve for N from this number.
assert(N==round(N));  % Make sure we got an integer.
if Pmin<0; Pmin=0; end
if Pmax>1; Pmax=1; end

%% Compute triangles scores
if ~isempty(scores_hist)
    cum_scores = [];
else
    [cum_scores, scores_hist, n_cores_1] = triangles_scores_mex...
        (N, Rij, hist_intervals, 0, 1, 0);

    % Normalize cumulated scores
    cum_scores = cum_scores/(2*(N-2));
    
    if verbose; log_message('Cores used for P computation: %d', n_cores_1); end
end

%% Histogram decomposition: P & sigma evaluation

% Initialize variables
h = 1/hist_intervals;
hist_x = ((h/2):h:(1-h/2))'; % x-values of the histogram
A = (N*(N-1)*(N-2)/2)/hist_intervals*(a+1); % normalization factor of one component of the histogram
start_values = [0, 0.5, 2.5, 0.78];%0.78]; % B, P, b, x0 % GGG try to guess initial x0 automatically by local minimum; try to fit b(x0) to remove DoF; consider normalize B explicitly to remove DoF.
B0 = start_values(2)*(N*(N-1)*(N-2)/2) /... % normalization of 2nd component: B = P*N_delta/sum(f), where f is the component formula
    sum(((1-hist_x).^start_values(3)).*exp(-start_values(3)/(1-start_values(4)).*(1-hist_x)));
start_values(1) = B0;

% Fit the distribution
conf = fitoptions('Method','NonlinearLeastSquares',...
    'Robust','LAR',...
    'Lower',[0, Pmin^3, 2, 0],... % B, P, b, x0
    'Upper',[Inf, Pmax^3, Inf, 1],... % GGG: [100*B0, Pmax^3, 1e2, 1],... % bounded domains to force convergence of fit()
    'Startpoint',start_values);
ft = fittype(['(1-P)*' num2str(A) '*(1-x)^' num2str(a) ' + P*B*((1-x)^b)*exp(-b/(1-x0)*(1-x))'],'options',conf);
[c,gof] = fit(hist_x, scores_hist, ft);

% Derive P and sigma
Rsquare = gof.rsquare;
P = c.P^(1/3);
peak = c.x0;
sigma = (1-peak)/peak2sigma;

if verbose
    log_message('Estimated CL Errors P,STD:\t%.1f%%\t%.1f\n\tR-squared:\t%.2f',100*P,sigma,Rsquare);
end

%% Local histograms analysis
% calculate the probability Pij of every relative rotation Rij to be
% indicative (rather than arbitrary)

% Initialize parameters for probabilities computations
A = a+1; % distribution 1st component normalization factor
B=c.B/((N*(N-1)*(N-2)/2)/hist_intervals); % distribution 2nd component normalization factor
b=c.b; x0=c.x0; % distribution parameters as empirically derived from the histogram

% Calculate probabilities
[ln_f_ind, ln_f_arb, n_cores_2] = pairs_probabilities_mex(N, Rij, P^2,A,a,B,b,x0, 0, 0);
Pij = 1 ./ (1 + (1-P)/P*exp(ln_f_arb-ln_f_ind));

% Fix singular output
if any(isnan(Pij))
    warning('NaN''s in probabilities computations uccured %d times out of %d!', sum(isnan(Pij)), numel(Pij));
    Pij(isnan(Pij)) = 0;
end

%% Summary
scores_fit.c = c;
scores_fit.gof = gof;
if verbose
    log_message('Cores used for Pij computation: %d', n_cores_2);
    log_message('Common lines probabilities to be indicative <Pij>=%.1f%%', 100*mean(Pij));
end

end
