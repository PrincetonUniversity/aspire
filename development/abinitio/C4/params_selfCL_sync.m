function params_simul = params_selfCL_sync

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params_simul.debug = true; % if true then analysis is performed against GT
params_simul.nImages = 100; % number of images
params_simul.c4_VolType = 'TRPV1'; % either GAUSSIAN,SYNTHETIC or IP3 or TRPV1
% params_simul.SNR = 1/16;
params_simul.SNR = 100000000;
% finding common-lines. Either 'GT' (ground-truth) or 'CALC' (calculate)
params_simul.CL  = 'CALC';
% finding self common-lines. Either 'GT' (ground-truth) or 'CALC' (calculate)
params_simul.SCL = 'CALC';

params_simul.EQUATOR = 'GT'; % can be either 'GT' or 'CALC'

% common-lines ambiguity simulation
if strcmp(params_simul.CL,'GT')
    % simulation of J ambiguity (handedness)
    params_simul.confuse_cl_J = true;
    % simulation of g ambiguity (group ambiguity)
    params_simul.confuse_cl_g = true;
end

% self common-lines ambiguity simulation
if strcmp(params_simul.SCL,'GT')
    % simulation of J ambiguity (handedness)
    params_simul.confuse_scl_J = true; % can be either true or false
    % simulation of first vs third self-common-line ambiguity
    params_simul.confuse_scl   = true;
end

end