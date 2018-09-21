if ~isoctave()
    mex -lgomp omp_get_max_threads_mx.c
    mex -lgomp omp_set_num_threads_mx.c
else
    mex -lgomp omp_get_max_threads_mx.c
    mex -lgomp omp_set_num_threads_mx.c
end
