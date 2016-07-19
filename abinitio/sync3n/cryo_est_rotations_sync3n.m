
function [V, data, algo, disp] = main_3N (conf)
% CryoEM: reconstruct a 3D volume by its 2D projections images
% Using the "3N X 3N" algorithm introduced in <GGG insert link or name>
% 
% Input:
% conf - configuration struct whose possible fields are summarized in
% initial_configuration.m
% 
% Output:
% V - reconstructed 3D volume
% data - struct containing variables describing input & output data
% algo - struct containing parameters of the algorithm
% disp - struct containing parameters of the display and saving of the results
% 
% Ido Greenberg, 2015

tic;

% Configuration
[data, algo, disp] = initial_configuration(conf);
if disp.SAVE_DATA
    save(sprintf('%s/configuration.mat', data.meta.PATH));
    diary(sprintf('%s/stdout.txt', data.meta.PATH));
end
fprintf('\n');
print_main_configurations(data);


% We wish to allow running the algorithm from the middle. This is
% particularly useful in two cases:
% 1. when we do several experiments with similar beginning, we can avoid
% unnecessary computations repetitions
% 2. when an experiment is surprisingly shut down
% 
% to support running from the middle, we check the type of the source file used
% for running, and skip parts of the algorithm accordingly.

if strcmp('SYNCHRONIZED_RIJ',data.meta.DATA_SOURCE)
    
    % keep necessary data before overwriting them by loaded data
    tmp_metadata = keep_necessary_data(data.meta);
    % load data from previous experiment
    load(data.meta.SOURCE_FILE);
    % restore the necessary data
    data.meta = restore_necessary_data(data.meta, tmp_metadata);

else
    
    if strcmp('UNSYNCHRONIZED_RIJ',data.meta.DATA_SOURCE)
    
        % keep necessary data before overwriting them by loaded data
        tmp_metadata = keep_necessary_data(data.meta);
        % load data from previous experiment
        load(data.meta.SOURCE_FILE);
        % restore the necessary data
        data.meta = restore_necessary_data(data.meta, tmp_metadata);
    
    else
        
        if strcmp('COMMON_LINES',data.meta.DATA_SOURCE)
            
            % keep necessary data before overwriting them by loaded data
            tmp_metadata = keep_necessary_data(data.meta);
            % load data from previous experiment
            load(data.meta.SOURCE_FILE);
            % restore the necessary data
            data.meta = restore_necessary_data(data.meta, tmp_metadata);
    
        else
            % Strat from the beginning
            
            % Create Data
            data = create_data (data, algo.N_THETA, disp.SAVE_DATA);
            print_time('Load Data',true);
            
            % Clear Images
            data.xxN1.projs = clear_images (...
                data.xxN1.projs_uncleared, algo.PRE_CLEAR, data.meta, disp.SAVE_DATA);
            
            data = save_and_cut(data, data.meta.N, 'Clear Images', disp.SAVE_DATA, []);
            
            % Common Lines
            [data.N1N1.clstack, data.N1N1.corrstack, ~, ~] = find_common_lines...
                (data.xxN1.projs, data.meta.n, algo.N_THETA, algo.shifts, false);
            % note that we don't use the shifts estimations returned here,
            % since we can better estimate them later (after estimating
            % rotations).
            
        end
        
        data = save_and_cut(data, algo.N_relative_rotations, 'Common Lines', disp.SAVE_DATA,...
            classificate_cls(data.N1N1.corrstack));

        % Estimate relative rotations
        [data.xxN2.Rij0, data.N2.r_valid_ks, data.N2.r_good_ks, data.N2.peak_width] = ...
            R_pairs(algo, data.N1N1.clstack, disp.R_pairs_progress_update);

    end
    
    data = save_and_cut(data, algo.N_triangles_sync, 'Relative Rotations', disp.SAVE_DATA,...
        classificate_Rijs(data.N1N1.corrstack, data.N2.r_valid_ks, data.N2.r_good_ks, data.N2.peak_width));
    
    % Reflections synchronization over relative rotations
    % local synchronization (for evaluation of relative rotations data reliability)
    data = triangles_sync(data, algo, disp.SAVE_DATA);
    print_time('Initial Triangles Synchronization',true);
    % global synchronization
    [data.N2.J_sync, data.N2.J_significance, data.others.J_sync_eigs, data.others.J_sync_iterations, data.others.J_sync_convergence] = ...
        signs_power_method(data.meta.N, data.xxN2.Rij0, algo.J_SYNC_N_EIGS, algo.SCORES_AS_ENTRIES, 2, false);
    % Do synchronize J-multiplications between relative rotations
    data.xxN2.Rij = sync_J(data.N2.J_sync, data.xxN2.Rij0);
    
end

% Save current data
save_and_cut(data, data.meta.N, 'Outer J Sync', disp.SAVE_DATA, []);
print_time('Outer J Synchronization',true);

% Build 3NX3N Matrix S
data.N3N.S = sync_matrix(data.meta.N, data.xxN2.Rij);

% Absolute Rotations
[data.xxN1.R, data.others.main_component, data.others.eigenvalues, data.N2.weights, data.others.sync_consistency] =...
    find_rotations(data.meta.N, algo, data.N3N.S, data.N2, data.N1N1, disp, data.meta.PATH);

data = save_and_cut(data, sum(data.others.main_component),...
    'Absolute Rotations Synchronization', disp.SAVE_DATA, data.others.main_component, false);


% Save
if disp.SAVE_DATA
    save_final_data (data, algo, disp, V);
end

% End
print_time('',true);
if disp.SAVE_DATA
    diary off;
end

end
