% OMP_SET_NUM_THREADS Set number of OpenMP threads
%
% Usage
%    omp_set_num_threads(n_threads);
%
% Input
%    n_threads: The maximum number of threads used by OpenMP.

function omp_set_num_threads(n_threads)
    if numel(n_threads) ~= 1 || ...
        abs(n_threads-round(n_threads)) > 1e-15 || ...
        n_threads < 1

        error('''n_threads'' must be a positive integer.');
    end

    try
        omp_set_num_threads_mx(n_threads);
    catch
        warning(['Cannot set number of OpenMP threads. ' ...
            'Try running ''install_mex''.']);
    end
end
