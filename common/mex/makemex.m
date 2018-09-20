if ~isoctave()
    mex CFLAGS="\$CFLAGS -fopenmp" -lgomp omp_get_num_threads_mx.c
    mex -lgomp omp_set_num_threads_mx.c
else
    mex omp_get_num_threads_mx.c
    mex omp_get_num_threads_mx.c
end
