function [stats,h]=cryo_analyze_sync_results(S,clstack,n_theta)
%
% CRYO_ANALYZE_SYNC_RESULTS   Analyze the results of the synchronization
%                             algorithm.
%
% cryo_analyze_sync_results(S,clstack,n_theta)
%       Compute statiscs for the rotations estimated by the syncronization
%       algorithm. The function takes the synchronization matrix S, the
%       common lines matrix clstack used to construct S, and n_theta, and
%       computes the statistics for the estimated rotations.
%       The statistics are computed by taking the estimated rotations
%       R_1,...,R_N computed from S and all common lines pairs c(i,j) and
%       c(j,i), 1<=j<i<=N and computing the angle R_i c(i,j) and R_j c(j,i)
%       for all i and j. These angles are the embedding errors of the
%       common lines. For this list of angles, we compute:
%           stats.mean      Mean embedding error
%           stats.std       STD of embedding error
%           stats.median    Median embedding error
%           stats.min       Min embedding error
%           stats.max       Max embedding error
%           stats.errs      List of embedding errors
%           stats.ev        Top 10 eigenvalues of S
%           stats.gap       ev(3)/ev(4).
%
% Yoel Shkolnisky, March 2015.

rotations=cryo_syncrotations(S);
ev=eigs(S,10);
ev=sort(ev,'descend');

clerr=cryo_syncconsistency(rotations,clstack,n_theta);
h=figure;
hist(clerr(:,3),360);
stats.mean=mean(clerr(:,3));
stats.std=std(clerr(:,3));
stats.median=median(clerr(:,3));
stats.min=min(clerr(:,3));
stats.max=max(clerr(:,3));
stats.err=clerr(:,3);
stats.ev=ev;
stats.gap=ev(3)/ev(4);