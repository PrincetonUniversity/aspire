
function [] = print_resolutions(n, Ns, NN, S_WEIGHTS, J_WEIGHTS, GROUPS)

% input
if ~exist('n','var') || isempty(n); n=89; end
if ~exist('Ns','var') || isempty(Ns); Ns=[1000]; end
if ~exist('NN','var') || isempty(NN); NN=[2]; end
if ~exist('S_WEIGHTS','var') || isempty(S_WEIGHTS); S_WEIGHTS=[1]; end
if ~exist('J_WEIGHTS','var') || isempty(J_WEIGHTS); J_WEIGHTS=[0]; end
if ~exist('GROUPS','var') || isempty(GROUPS); GROUPS=[1]; end

% enviromental params
BASE_PATH = '/home/idog/matlab/aspire_merge/output';
DATA_SET = sprintf('80s_%d',n);

% iterate
fprintf('\nN\tnn\tS_w\tJ_w\tg\tresolution [A]\n\n');
for N = Ns
    for nn = NN
        for S_weights = S_WEIGHTS
            for J_weights = J_WEIGHTS
                for g = GROUPS
                    serial = sprintf('N%d_nn%02d_sw%d_jw%d_g%d',...
                        N, nn, S_weights, J_weights, g);
                    res = load(sprintf('%s/%s/%s/resolution.mat',...
                        BASE_PATH, DATA_SET, serial));
                    fprintf('N=%d\tnn%02d\tsw=%d\tjw=%d\tg%d\t%.2f\n',...
                        N, nn, S_weights, J_weights, g, res.res);
                end
            end
        end
        fprintf('\n');
    end
end

end
