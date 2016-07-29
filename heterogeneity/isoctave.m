% ISOCTAVE Determine whether running on Octave
%
% Usage
%    ret = isoctave();
%
% Output
%    ret: True if system is Octave, false if MATLAB.

function ret = isoctave()
    persistent p_ret;
    if isempty(p_ret)
        p_ret = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
    end
    ret = p_ret;
end
