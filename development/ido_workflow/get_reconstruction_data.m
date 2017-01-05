
function data = get_reconstruction_data(n, data_set, N, nn, S_weights, J_weights, g)

if ~exist('n','var') || isempty(n); n=89; end
if ~exist('data_set','var') || isempty(data_set); data_set=sprintf('80s_%d',n); end
if ~exist('N','var') || isempty(N); N=1000; end
if ~exist('nn','var') || isempty(nn); nn=2; end
if ~exist('S_weights','var') || isempty(S_weights); S_weights=1; end
if ~exist('J_weights','var') || isempty(J_weights); J_weights=0; end
if ~exist('g','var') || isempty(g); g=1; end

BASE_PATH = '/home/idog/matlab/aspire_merge/output';
serial = sprintf('N%d_nn%02d_sw%d_jw%d_g%d',...
    N, nn, S_weights, J_weights, g);
base_path = sprintf('%s/%s/%s', BASE_PATH, data_set, serial);

data.conf = load(sprintf('%s/configuration', base_path));
data.projs = load(sprintf('%s/projections', base_path));
data.cls = load(sprintf('%s/common_lines', base_path));
data.rij = load(sprintf('%s/relative_rotations', base_path));
data.jsync = load(sprintf('%s/J_sync', base_path));
data.w = load(sprintf('%s/weights', base_path));
data.ri = load(sprintf('%s/absolute_rotations', base_path));
data.shifts = load(sprintf('%s/shifts_and_phaseflip', base_path));
data.fsc = load(sprintf('%s/resolution', base_path));

end
