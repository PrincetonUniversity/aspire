% OMP_GET_NUM_THREADS Get number of OpenMP threads
%
% Usage
%    n_threads = omp_get_num_threads();
%
% Output
%    n_threads: The maximum number of threads used by OpenMP.

function n_threads = omp_get_num_threads()
    try
        n_threads = omp_get_num_threads_mx();
    catch
        warning(['Cannot get number of OpenMP threads. ' ...
            'Try running ''install_mex''.']);
    end
end
