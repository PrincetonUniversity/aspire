function [J_sync,J_significance,eigenvalues,itr,dd] =...
    signs_power_method ( N , J_list , n_eigs , verbose , measure_time,s )
    %SIGNS_POWER_METHOD Calculate first eigenvalues & eigenvectors of the
    % Signs matrix using power method, the D2 version. 
    % 
    % Every iteration, the multiplication Signs*v is computed within an
    % internal function, which calculates Signs entries repeatedly by
    % triangles synchronization (to avoid memory over-usage).
    % 
    % See signs_times_v_mex() for more details.
    % 
    % Input:
    % N = number of projections on which Rij are based
    % Rij = 3X3XN-choose-2 array of estimations of {Ri'Rj}, before J-synchronization
    % n_eigs = number of eigenvalues to compute
    % score_as_entries = whether to use triplets scores as the matrix
    %  entries or just the signs
    % verbose = 0 for nothing, 1 for global printing, 2 for iterative printing
    % 
    % Output:
    % J_sync = N-choose-2 dimensional vector. each entry idx expresses the
    %  J-synchronization of the corresponding relative rotation Rij(:,:,idx)
    % J_significance = N-choose-2 dimensional vector. each entry idx expresses the
    %  significance of the corresponding relative rotation J-synchronization,
    %  based on the eigenvector of the signs matrix
    % eigenvalues = all the eigenvalues of the signs matrix that were computed
    % itr = number of iterations needed for the power method
    % dd = size of last change in the iterative estimated eigenvector
    % 
    % Written by Ido Greenberg, 2015
    % Adapted for D2 symmetry by Eitan Rosen, 2017. 
   % signs_power_method(nImages, reshape(qi_qjs,[3 3 nchoosek(nImages,2)]),1,0);
    % input validation
    if (n_eigs <= 0); error('n_eigs must be positive integer!'); end
    if ~exist('verbose','var'); verbose=0; end
    if ~exist('measure_time','var'); measure_time=false; end
    
    % constants
    epsilon = 0.0000;
    tol=1e-4;
    MAX_ITERATIONS = 100;
    if verbose >= 2
        fprintf('Power method for signs matrix started, up to eigenvector\naccuracy goal of %f, and limited by %d iterations.\n',...
            epsilon, MAX_ITERATIONS);
    end
    
    % initialization
    if exist('s','var')
        rng(s);
    end
    N_pairs = nchoosek(N,2);
    vec=rand(N_pairs,n_eigs);
    vec=vec/norm(vec);
    dd = 1;
    eig_diff=inf;
    itr = 0;
    eigenvalues=zeros(3,3);
    prev_eval=inf;
    
    % power method iterations, stopping condition also includes proximity
    % to the theoretic largest eigenvalue 2N-4. 
    if measure_time; t = toc; end
    while itr < MAX_ITERATIONS && ...
            (abs(eigenvalues(1,1)/(2*N-4))>1+epsilon || ...
                abs(eigenvalues(1,1)/(2*N-4))<1-epsilon) && ...
                    eig_diff>tol
                
        itr = itr + 1;
        %fprintf('itr no. %d\n',itr);
        vec_new = signs_times_v_mex2(N,J_list,vec,n_eigs);
        vec_new = reshape(vec_new,size(vec));
        [vec_new,eigenvalues] = qr(vec_new,0);
        dd = norm(vec_new(:,1)-vec(:,1));
        vec = vec_new;
        eig_diff=abs(prev_eval-eigenvalues(1,1));
        prev_eval=eigenvalues(1,1);
        if verbose >= 2
            fprintf('Iteration %02d: ||curr_evec-last_evec|| = %.3f\n', itr, dd);
        end
    end
    if measure_time; t=(toc-t)/60; end
    disp(['num of Jsync iterations: ' num2str(itr)]);
    % set output format
    eigenvalues = diag(eigenvalues);
    J_significance = abs(vec(:,1));
    J_sync = sign(vec(:,1)); % only signs of the eigenvector entries are needed
    
    % print eigenvalues & number of iterations needed for the power method
    if verbose >= 1
        fprintf('Outer J Synchronization:\n\titerations needed: %d\n', itr);
        if measure_time; fprintf('\ttime consumed: %f [min]\n', t); end
        fprintf('\tsigns matrix eigenvalues (in practice vs. theoretical):\n\t\t%.0f\t(%.0f)\n', eigenvalues(1), 2*N-4);
        for i=2:min(n_eigs,5); fprintf('\t\t%.0f\t(%.0f)\n', eigenvalues(i), N-4); end
    end

end
