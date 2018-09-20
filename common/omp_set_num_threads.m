% OMP_SET_NUM_THREADS Set number of OpenMP threads
%
% Usage
%    omp_set_num_threads(n_threads);
%
% Input
%    n_threads: The maximum number of threads used by OpenMP.

function omp_set_num_threads(n_threads)
    try
        omp_set_num_threads_mx(n_threads);
    catch
        warning(['Cannot set number of OpenMP threads. ' ...
            'Try running ''install_mex''.']);
    end
end
