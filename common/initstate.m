function initstate(state)
%
% Initialize the state of the random numbers generator.
% This enables reproducing the generated "random" numbers in different
% calls to rand and randn.
%
% Yoel Shkolnisky, January 2007.
% Revised: Yoel Shkolnisky, September 2013

if nargin<1
    state = 3498;
end

if ~isoctave()
    rng('default');
    rng(state);
else
    % There doesn't seem to be an overall way of resetting all the RNGs in
    % Octave at the same time so let's do it manually.

    rand('state', state);
    randn('state', state);
    rande('state', state);
    randg('state', state);
    randp('state', state);
end

