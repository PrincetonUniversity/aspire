function cryo_workflow_abinitio(workflow_fname)
% CRYO_WORKFLOW_ABINITIO  Interactive generation of abinitio models
%
% cryo_workflow_classify(workflow_fname)
%   Interactively collect all parameters required to generate abinitio
%   models from precomputed class averages and execute the reconstruction.
%   workflow_name is the name of the file where the entered parameters will
%   be saved.
%
% See also cryo_workflow_start
%
% Yoel Shkolnisky, August 2015.

% Read workflow file
tree=xmltree(workflow_fname);
workflow=convert(tree);

%% Read Classification parameters

defnmeans=1000; % Default number of means to use to abinitio reconstruction.
message=sprintf('Number of class means to use to abinitio reconstruction of each group? ');
nmeans=fmtinput(message,defnmeans,'%d');

% Find all available class means files.
avgfiles=dir(fullfile(workflow.info.working_dir,'averages_nn*_group*.mrc'));

available_nns=zeros(numel(avgfiles),1);
for k=1:numel(avgfiles)
    avg_file_params=sscanf(avgfiles(k).name,'averages_nn%d_group%d.mrc');
    available_nns(k)=avg_file_params(1);
end
available_nns=sort(unique(available_nns)); % Available nnavg values.

% Convert available_nns into a cell array for use as input in
% multichoice_question.
cell_nn=cell(numel(available_nns),1);
for k=1:numel(available_nns)
    cell_nn{k}=num2str(available_nns(k));
end
message=sprintf('Averaging value (nnavg) to use ');
nnavg=multichoice_question(message,cell_nn,available_nns,num2str(available_nns(end)));

% Choose reconstruction algorithm
message='Which abinitio reconstruction algorithm to use?';
algo=multichoice_question(message,{'sync3N','sync2N','LUD'},[ 1, 2, 3],'sync3N');


%% Update workflow struct

workflow.abinitio.nmeans=nmeans; % Number of means to use to abinitio reconstruction.
workflow.abinitio.nnavg=nnavg;   % nnavg to use to abinitio reconstruction.
workflow.abinitio.algo=algo;

tree=struct2xml(workflow);
save(tree,workflow_fname); 

%% Generate abinitio models
cryo_workflow_abinitio_execute(workflow_fname);
