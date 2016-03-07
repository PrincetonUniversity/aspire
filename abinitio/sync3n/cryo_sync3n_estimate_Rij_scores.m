function [P, sigma, Rsquare, Pij, scores_hist, cum_scores, plotted] =...
    cryo_sync3n_estimate_Rij_scores(Rij,required_P, PLOT, verbose)
% CRYO_SYNC3N_ESTIMATE_RIJ_SCORES   Reliability score of relative rotations.
%
% [P, sigma, Rsquare, Pij, scores_hist, cum_scores, plotted] =...
%   cryo_sync3n_estimate_Rij_scores(Rij,required_P, PLOT, verbose)
%       Given a list of relative rotations between pairs of images,
%       computed for example using cryo_sync3n_estimate_all_Rijs, estimate
%       several reliability indicators for each relative rotations in the
%       list.
% 
% Input:
% Rij           Relative rotations matrices (3 X 3 X N-choose-2)
% required_P    Minimum required value for estimated rate of indicative
%       common lines.  If rate is estimated to be lower than required_P,
%       the function aborts (to save time) and the caller should improve
%       the data and call again. Set to 0 to allow any P. Default value is
%       XXX.
% PLOT  Whether to plot the histogram & its fit. Default XXX. 
%       XXX What is actually plotted?
% verbose Set to nonzero to print summary messages. Default 0.
% 
% Output:
% P         Estimated rate of indicative common lines. 
%       XXX Somewhere you should define indicative.
% sigma     Typical error in indicative common lines estimations [degrees].
% Rsquare   R^2 value of the fit of the histogram which derives P & sigma
% Pij       Estimated probabilities of each common lines to be indicative
% scores_hist   Histogram of the common lines triangles scores 
%           XXX What is this (figure handle? array?)
% cum_scores    Pairs of cumulative triangles scores
% plotted   Whether the histogram was plotted or not (depends on both PLOT
%           and the flow of the function)
%           XXX Why do you need this? Can't you use the returned figure
%           handle to indicate that?
% 
% This is a wrapper for the mex functions triangles_scores_mex and
% pairs_probabilities_mex.
% More information can be found in the mex functions.
%
% This function was revised from 
%   triangles_scores(Rij, N, required_P, PLOT, verbose)
% Written by Ido Greenberg, 2015
% Revised by Yoel Shkolnisky, March 2016.

if ~exist('verbose','var')
    verbose=0;
end

N=(1+sqrt(1+8*size(Rij,3)))/2; % Extract N from the number of relative rotaitons in Rij.
% The number of relative rotation for a given N is (N choose 2) so,
% solve for N from this number.
assert(N==round(N));  % Make sure we got an integer.

% Histogram Analysis Constants
hist_intervals = 100; % number of intervals in the scores histogram
a = 2.0; % empirical constant factor in the theoretical formula of the scores distribution
% GGG 'a' was 2.0 in the original version, but 2.5 seems better, at least for decomposition_test with N~100
peak2sigma = 2.43e-2; % empirical relation between the location of the peak of the histogram, and the mean error in the common lines estimations

%% Compute triangles scores

[cum_scores, scores_hist] = triangles_scores_mex(N, Rij, hist_intervals, 0);

% Normalize cumulated scores
%cum_scores = ( cum_scores/(2*(N-2)) ).^0.5; % mean squares
cum_scores = cum_scores/(2*(N-2)); % mean

%% Histogram decomposition: P & sigma evaluation

% Initialize variables
h = 1/hist_intervals;
hist_x = ((h/2):h:(1-h/2))'; % x-values of the histogram
A = (N*(N-1)*(N-2)/2)/hist_intervals*(a+1); % normalization factor of one component of the histogram
start_values = [0, 0.5, 2.5, 0.9]; % B, P, b, x0
B0 = start_values(2)*(N*(N-1)*(N-2)/2) /... % normalization of 2nd component: B = P*N_delta/sum(f), where f is the component formula
    sum(((1-hist_x).^start_values(3)).*exp(-start_values(3)/(1-start_values(4)).*(1-hist_x)));
start_values(1) = B0;

% Fit the distribution
conf = fitoptions('Method','NonlinearLeastSquares',...
    'Robust','LAR',...
    'Lower',[0, 0, 2, 0],... % B, P, b, x0
    'Upper',[Inf, 1, Inf, 1],...
    'Startpoint',start_values);
ft = fittype(['(1-P)*' num2str(A) '*(1-x)^' num2str(a) ' + P*B*((1-x)^b)*exp(-b/(1-x0)*(1-x))'],'options',conf);
[c,gof] = fit(hist_x, scores_hist, ft);

% Derive P and sigma
Rsquare = gof.rsquare;
P = c.P^(1/3);
peak = c.x0;
sigma = (1-peak)/peak2sigma;

% Check validity of estimated P, and plot histogram if needed
plotted = false;
if PLOT && P >= required_P
    figure;
    title('Triangles Scores');
    hold on;
    h = bar(hist_x,scores_hist,'facecolor','c');
    %set(get(h,'child'),'facea',0.4);
    hplot = plot(c);
    set(hplot, 'linewidth',1.5);
    xlabel('');
    ylabel('');
    legend('hide');
    plotted = true;
end

%% Local histograms analysis
% calculate the probability Pij of every relative rotation Rij to be
% indicative (rather than arbitrary)

% Initialize parameters for probabilities computations
A = a+1; % distribution 1st component normalization factor
B=c.B/((N*(N-1)*(N-2)/2)/hist_intervals); % distribution 2nd component normalization factor
b=c.b; x0=c.x0; % distribution parameters as empirically derived from the histogram

% Calculate probabilities
[ln_f_ind, ln_f_arb] = pairs_probabilities_mex(N, Rij, 0, P^2,A,a,B,b,x0);
Pij = 1 ./ (1 + (1-P)/P*exp(ln_f_arb-ln_f_ind));

% Fix singular output
if any(isnan(Pij))
    warning('NaN''s in probabilities computations uccured %d times out of %d!', sum(isnan(Pij)), numel(Pij));
    Pij(isnan(Pij)) = 0;
end

%% Summary
if verbose
    log_message('Estimated CL Errors P,STD:\t%.1f%%\t%.1f\n\tR-squared:\t%.2f',100*P,sigma,Rsquare);
    log_message('Percentage of common lines with Pij>0.5: %.1f%% (Pij is the probablity of the common lines between i and j to be indicative).', 100*sum(Pij>0.5)/numel(Pij));
end

end
