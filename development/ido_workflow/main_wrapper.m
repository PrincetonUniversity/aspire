
function [] = main_wrapper(n, DATA_SET, Ns, NN, S_WEIGHTS, J_WEIGHTS, GROUPS)

% input
if ~exist('n','var') || isempty(n); n=89; end
if ~exist('DATA_SET','var') || isempty(DATA_SET); DATA_SET=sprintf('80s_%d',n); end
if ~exist('Ns','var') || isempty(Ns); Ns=[1000]; end
if ~exist('NN','var') || isempty(NN); NN=[2]; end
if ~exist('S_WEIGHTS','var') || isempty(S_WEIGHTS); S_WEIGHTS=[1]; end
if ~exist('J_WEIGHTS','var') || isempty(J_WEIGHTS); J_WEIGHTS=[1]; end
if ~exist('GROUPS','var') || isempty(GROUPS); GROUPS=[1]; end

% environmental params
BASE_PATH = '/home/idog/matlab/results/cls_map_filter';

% general conf
MOLEC_RADIUS = 0.3;
N_THETA = 144;
FILTER_RADIUS = 4;
MAX_SHIFT = round(6*n/89);
SHIFT_STEP = max(ceil(MAX_SHIFT/20),1);
VOTING_TICS_WIDTH = 1;
P_PERMITTED_INCONSISTENCY = 1.5;
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
                    base_path = sprintf('%s/%s/%s', BASE_PATH, DATA_SET, serial);
                    mkdir(base_path);
                    
                    % run
                    main_projs2v (...
                        source_file, base_path, N, MOLEC_RADIUS,...
                        N_THETA, FILTER_RADIUS, MAX_SHIFT, SHIFT_STEP,...
                        VOTING_TICS_WIDTH, J_weights, S_weights,...
                        P_PERMITTED_INCONSISTENCY, REF_VOL);
                    
                end
            end
        end
    end
end
