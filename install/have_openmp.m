% HAVE_OPENMP Determines whether OpenMP is present on the system
%
% Usage
%    have_library = have_openmp();
%
% Output
%    have_library: True if OpenMP is installed and configured for either the
%       'gcc' or 'icc' compilers.
%
% Description
%   Checks whether 'gcc' or 'icc' are able to compile a simple OpenMP C
%   program and tries to run it. If successful, returns true. Note that
%   OpenMP may still be available on the system through other compilers.

function have_library = have_openmp()
    compilers = {'gcc', 'icc'};
    omp_flags = {'-fopenmp', '-qopenmp'};

    test_dir = tempname();

    mkdir(test_dir);

    test_c = fullfile(test_dir, 'test.c');
    test_exe = fullfile(test_dir, 'test');

    fd = fopen(test_c, 'w');
    fprintf(fd, '%s', [ ...
        '#include <stdio.h>' char(10) ...
        'char omp_set_num_threads ();' char(10) ...
        'int main ()' char(10) ...
        '{' char(10) ...
        'printf("%hhu\n", omp_get_num_threads ());' char(10) ...
        'return 0;' char(10) ...
        '}' char(10)]);
    fclose(fd);

    have_library = false;

    for k = 1:numel(compilers)
        cmd = sprintf('%s --version', compilers{k});
        [ret, ~] = system([cmd ' 2>&1']);

        if ret ~= 0
            continue;
        end

        cmd = sprintf('%s %s %s -o %s', ...
            compilers{k}, omp_flags{k}, test_c, test_exe);
        [ret, val] = system([cmd ' 2>&1']);

        if ret ~= 0
            continue;
        end

        cmd = test_exe;
        [ret, ~] = system([cmd ' 2>&1']);

        have_library = have_library | (ret == 0);

        delete(test_exe);
    end

    delete(test_c);
    rmdir(test_dir);
end
