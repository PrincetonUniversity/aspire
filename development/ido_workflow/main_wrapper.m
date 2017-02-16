
function [] = main_wrapper(n, ALGO, Ns, NN, S_WEIGHTS, J_WEIGHTS, GROUPS)

% input
if ~exist('n','var') || isempty(n); n=89; end
if ~exist('DATA_SET','var') || isempty(ALGO); ALGO='weighting_paper'; end
if ~exist('Ns','var') || isempty(Ns); Ns=[1000]; end
if ~exist('NN','var') || isempty(NN); NN=[2]; end
if ~exist('S_WEIGHTS','var') || isempty(S_WEIGHTS); S_WEIGHTS=[1]; end
if ~exist('J_WEIGHTS','var') || isempty(J_WEIGHTS); J_WEIGHTS=[0]; end
if ~exist('GROUPS','var') || isempty(GROUPS); GROUPS=[1]; end

% environmental params
%BASE_PATH = '/home/idog/matlab/aspire_merge/output';
BASE_PATH = sprintf('/home/idog/matlab/results/%s/80s_%d',ALGO,n);

% general conf
N_THETA = 360;
MOLEC_RADIUS = 0.4;
SHIFT_STEP = 1;
MAX_SHIFT = floor(10*n/89);
VOTING_TICS_WIDTH = 1;
REF_VOL = sprintf('/home/idog/matlab/aspire_merge/output/vol_80s_%d.mrc', n);

% iterate
for N = Ns
    for nn = NN
        for S_weights = S_WEIGHTS
            for J_weights = J_WEIGHTS
                for g = GROUPS
                
                    close all;
                    
                    % conf
                    source_file = sprintf(...
                        '/home/idog/matlab/dbs/80s_%d/averages_nn%02d_group%d.mrc',...
                        n, nn, g);
                    serial = sprintf('N%d_nn%02d_sw%d_jw%d_g%d',...
                        N, nn, S_weights, J_weights, g);
                    base_path = sprintf('%s/%s', BASE_PATH, serial);
                    mkdir(base_path);
                    
                    % run
                    main_projs2v (...
                        source_file, base_path, N,...
                        MOLEC_RADIUS, N_THETA, MAX_SHIFT, SHIFT_STEP,...
                        VOTING_TICS_WIDTH, [], J_weights, S_weights,...
                        REF_VOL);
                    
                end
            end
        end
    end
end
