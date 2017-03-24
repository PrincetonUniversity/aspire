
% tester

%main_wrapper(n, DATA_SET, Ns, NN, S_WEIGHTS, J_WEIGHTS, GROUPS)

% different confs
main_wrapper(89, [], 200, [7,5], [0,1], [0,1], []);
% nn02 paper conf
main_wrapper(89, [], 1000, [2], [1], [0], []);

% print results
print_resolutions(89, 200, [7,5], [0,1], [0,1], []);
print_resolutions(89, 1000, [2], [1], [0], []);
