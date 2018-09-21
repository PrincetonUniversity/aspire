if ~isoctave()
    mex CFLAGS="\$CFLAGS -fopenmp" -lgomp omp_get_max_threads_mx.c
    mex -lgomp omp_set_num_threads_mx.c
else
    env_cflags = getenv('CFLAGS');

    oct_cflags = mkoctfile('-p', 'CFLAGS');
    oct_cflags = [oct_cflags(1:end-1) ' -fopenmp'];

    setenv('CFLAGS', oct_cflags);

    mex -lgomp omp_get_max_threads_mx.c
    mex -lgomp omp_set_num_threads_mx.c

    setenv('CFLAGS', env_cflags);
end
