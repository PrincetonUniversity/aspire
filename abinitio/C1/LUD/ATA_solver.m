function A=ATA_solver(V1,V2,K)
% We look for a linear transformation (3 x 3 matrix) A such that
    % V1*A'=R1 and V2*A'=R2 are the columns of the rotations matrices
    %
    % Therefore:
    %
    % V1 * A'*A V1' = 1
    % V2 * A'*A V2' = 1
    % V1 * A'*A V2' = 0
    %
    % These are 3*K linear equations for to 9 matrix entries of A'*A
    % Actually, there are only 6 unknown variables, because A'*A is symmetric
    
    % 3*K equations in 9 variables (3 x 3 matrix entries)
    equations = zeros(3*K,9);
    
    for i=1:3;
        for j=1:3;
            equations(1:3:3*K, 3*(i-1)+j) = V1(i,:) .* V1(j,:);
            equations(2:3:3*K, 3*(i-1)+j) = V2(i,:) .* V2(j,:);
            equations(3:3:3*K, 3*(i-1)+j) = V1(i,:) .* V2(j,:);
        end;
    end;
    
    % Truncate from 9 variables to 6 variables corresponding
    % to the upper half of the matrix A'*A
    %truncated_equations = equations(:, [1, 2, 3, 5, 6, 9]);
    
    % b = [1 1 0 1 1 0 ...]' is the right hand side vector
    b = ones(3*K,1);
    b(3:3:3*K) = zeros(K,1);
%     alpha=sqrt(0.95*K);
    
    % Find the least squares approximation
    ATA = equations\b;

    ATA = reshape(ATA, 3*ones(1, 2));
   
    if any(eig(ATA)<=0)
        error('ATA is not positive definite');
    end

    % The Cholesky decomposition of A'*A gives A
    A = chol(ATA);
