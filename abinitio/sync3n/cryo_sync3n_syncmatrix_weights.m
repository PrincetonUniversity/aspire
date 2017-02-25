
function [W, Pij, scores_hist, cum_scores] = cryo_sync3n_syncmatrix_weights (Rij)
%TRIANGLES_SYNC synchronize J-multiplications between all triplets of
% relative rotations {Rij,Rjk,Rki}, and evaluate the probability of every
% relative rotation to be indicative ("good", rather than arbitrary).
% The probabilities are formed as a matrix of weights, used later for weighting
% the blocks {Rij} of the synchronization matrix S.
% 
% $$$ A reference to the weighting paper should appear here.
% 
% The probabilities are calculated in triangles_scores() (see documentation
% inside). This function is an iterative wrapper for these calculations.
% 
% - Why do the calculations have to be repeated iteratively?
% triangles_scores() first estimates the rate P of indicative relative
% rotations. Then it uses Bayesian approach to find the probability Pij of
% every relative rotation Rij to be indicative.
% When P is correctly estimated, then it is expected to approximately equal
% mean(Pij). Unfortunately, when P<<1, the estimation of P has a strong
% bias upwards. In this case {Pij} turn out to be far too small, and we
% have mean(Pij)<<P. When we detect this phenomenon, we try to re-estimate
% P, forcing stricter limits on the estimation. We do that iteratively,
% until we have P~mean(Pij).
% 
% Input:
%   Rij (3x3x(N-choose-2)): relative rotations between viewing directions.
% 
% Output:
%   W (NxN): weights matrix (later we will have S <- S.*W).
%   Pij (N-choose-2): the probabilities of the relative rotations to be
%     indicative.
%   scores_hist (struct): properties of the histogram of the triangles scores
%     (used for the estimation of Pij). P is included here, and can be used
%     as indication for the quality of the data and the reliability of the
%     reconstruction.
%   cum_scores (N-choose-2): cummulated triangles scores - sum of scores
%     {s_ijk} for every pair ij.
% 
% Written by Ido Greenberg, 2016

%% Configuration
PERMITTED_INCONSISTENCY = 1.5; % Consistency condition is: mean(Pij)/PERMITTED_INCONSISTENCY < P < mean(Pij)*PERMITTED_INCONSISTENCY.
P_DOMAIN_LIMIT = 0.7; % Forced domain of P is [Pmin,Pmax], with Pmin=P_DOMAIN_LIMIT*Pmax.
MAX_ITERATIONS = 12; % Maximum iterations for P estimation.
MIN_P_PERMITTED = 0.04; % When P is that small, stop trying to synchronize P with Pij, since we have no chance.

%% Get initial estimation for Pij

% get estimation
scores_hist = struct;
[P, scores_hist.sigma, scores_hist.Rsquare, Pij, scores_hist.hist, scores_hist.fit, cum_scores] =...
    triangles_scores(Rij);

% are P and Pij consistent?
too_low = P < mean(Pij)/PERMITTED_INCONSISTENCY;
too_high = P > mean(Pij)*PERMITTED_INCONSISTENCY;
inconsistent = too_low | too_high;

% define limits for next P estimation
if too_high
    Pmax = P;
    Pmin = Pmax*P_DOMAIN_LIMIT;
else
    Pmin = P;
    Pmax = Pmin/P_DOMAIN_LIMIT;
end

%% Repeat iteratively until estimations of P & Pij are consistent

iteration = 0;

while inconsistent
    
    % count iterations
    iteration = iteration + 1;
    if iteration >= MAX_ITERATIONS
        warning('Triangles Scores are too bad distributed, P failed to converge wrt Pij.');
        break;
    end
    
    % keep last inconsistency direction to recognize trends
    prev_too_low = too_low;
    %prev_too_high = too_high;
    
    % get estimation
    [P, scores_hist.sigma, scores_hist.Rsquare, Pij, ~, scores_hist.fit, ~] = ...
        triangles_scores(Rij, scores_hist.hist, Pmin, Pmax);
    
    % are P and Pij consistent?
    too_low = P < mean(Pij)/PERMITTED_INCONSISTENCY;
    too_high = P > mean(Pij)*PERMITTED_INCONSISTENCY;
    inconsistent = too_low | too_high;
    
    % check trend
    if too_low ~= prev_too_low
        % error trend inverted -> we're around the correct value -> use stricter limits
        P_DOMAIN_LIMIT = sqrt(P_DOMAIN_LIMIT);
    end
    
    % define limits for next P estimation
    if too_high
        if P < MIN_P_PERMITTED
            error('Triangles Scores are too bad distributed, whatever small P we force.');
        end
        Pmax = Pmax*P_DOMAIN_LIMIT;
        Pmin = Pmax*P_DOMAIN_LIMIT;
    else
        Pmin = Pmin/P_DOMAIN_LIMIT;
        Pmax = Pmin/P_DOMAIN_LIMIT;
    end
    
end

scores_hist.P = P;

%% Fill weights matrix

N = 0.5*(1+sqrt(1+8*size(Rij,3)));
W = zeros(N);
idx = 0; % pair index

for i = 1:N
    for j = (i+1):N
        idx = idx + 1;
        W(i,j) = Pij(idx);
        W(j,i) = Pij(idx);
    end
end

end
