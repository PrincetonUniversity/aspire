function sign_J = cryo_outer_sync_Rij(Rijs, K)

%
% Synchronize the J-duplicity on the graph consisting of nodes of all
% (i,j)'s.
%
% input: Rijs, the estimated R_i^TR_j up to J-duplicity, n_node = K*(K-1)/2
%
% output:   Rotations, 3-by-3-by-K, the R_i. The 3-top eigenvecotors of bigS
%           upto a normalization constant are [R_1^T; R_2^T;...; R_K^T]
%
% Y. Shkolnisky and X. Cheng, July 2012.

if (size(Rijs,1)~=9) || (size(Rijs,2)~=K*(K-1)/2)
    error('Rijs must be 9xK(K-1)/2');
end

n = K;

%% build the J-synchronization matrix S

n_edge = n*(n-1)*(n-2)/2; % the number of triangles is C(3,n). Every triangle induces 3 edges.
II = zeros(n_edge,1);
JJ = zeros(n_edge,1);
VV = zeros(n_edge,1);

poolreopen;

for i=1:n
    % j spans from (i+1) to (n-1). Hence n-i-2 elements are reserved
    II_temp = cell(1,n-i-2);
    JJ_temp = cell(1,n-i-2);
    VV_temp = cell(1,n-i-2);
    
    parfor j=(i+1):(n-1)
        [II_temp{j},JJ_temp{j},VV_temp{j}] = computeTriangles(Rijs,i,j,n)
    end
    
    % jump directly to the next entry : the number of entries already induced 
    % by i'=1,2,...,i-1 is 3*(C(2,n-1)+C(2,n-2)+...+C(2,n-(i-1))) = 3*(C(3,n)-C(3,n-i+1))
    ii=(n*(n-1)*(n-2)-(n-i+1)*(n-i)*(n-i-1))/2;
    
    for j=(i+1):(n-1)
        %every j induces 3*(n-j) k's
        inds = ii+1:ii+3*(n-j);
        II(inds) = II_temp{j};
        JJ(inds) = JJ_temp{j};
        VV(inds) = VV_temp{j};
        ii = ii + 3*(n-j);
    end
end

n_node = n*(n-1)/2;
S = sparse(II,JJ,VV, n_node, n_node);
S = S+ S';

% since we are only interested in the first eigenvector, power method is
% more efficient
[~,evect1,~] = powerMethod(S);


%  figure(1),
%  subplot(2,1,1); hist(evals, 40);
%  set(gca, 'FontSize', 20);
%  title('Outer Sync');



% check spectrum...


%% build bigS and reconstruct the rotations

% make J-conjugate on Rij's indicating by the top eigenvector of S
sign_J = sign(evect1); %indicating if R_(i,j) should be J-conjugated;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [II,JJ,VV] = computeTriangles(Rijs,i,j,n)


% if Rij=Ri^TRj for i,j = 1,2,3, then the three equations
%   Rij*Rjk=Rki.',    Rjk*Rki=Rij.',    Rki*Rij=Rjk.'
% are satisfied, and the function returns 0. Otherwise it returns a
% positive value which is the "discrepancy".
% discrepancyRij=@ (R12, R23, R31)...
%     norm(R12*R23 - R31.', 'fro')^2 +norm(R23*R31 - R12.', 'fro')^2+norm(R31*R12 - R23.', 'fro')^2;

J = diag([1 1 -1]); %flipping matrix

n_edge = 3*(n-j);
II = zeros(n_edge,1);
JJ = zeros(n_edge,1);
VV = zeros(n_edge,1);

ind1 = uppertri_ijtoind(i,j,n);
R12 = reshape(Rijs(:,ind1), 3,3);
ii = 1;

for k=j+1:n
    ind2 = uppertri_ijtoind(i,k,n);
    ind3 = uppertri_ijtoind(j,k,n);
    R23 = reshape(Rijs(:,ind3), 3,3);
    R31 = reshape(Rijs(:,ind2), 3,3).';
    
    % compute 4 discrepancyRij values, the lowest one indicating the
    % +-1 on the edges of the triangle (ij)-(jk)-(ki)
    y = zeros(1,4);
    y(1) = discrepancyRij(R12, R23, R31);
    y(2) = discrepancyRij(J*R12*J, R23, R31);
    y(3) = discrepancyRij(R12, J*R23*J, R31);
    y(4) = discrepancyRij(R12, R23, J*R31*J);
    
    [~, imin] = min(y);
    if imin == 1
        vals_edge = [1 1 1];
    else if imin == 2
            vals_edge = [-1 1 -1];
        else if imin == 3
                vals_edge = [-1 -1 1];
            else
                vals_edge = [1 -1 -1];
            end
        end
    end
    
    II(ii) = ind1;  JJ(ii) = ind3;  VV(ii)=vals_edge(1); ii=ii+1;
    II(ii) = ind3;  JJ(ii) = ind2;  VV(ii)=vals_edge(2); ii=ii+1;
    II(ii) = ind2;  JJ(ii) = ind1;  VV(ii)=vals_edge(3); ii=ii+1;
end

end


function discrepancy = discrepancyRij(R12, R23, R31)
    discrepancy = ...
        norm(R12*R23 - R31.', 'fro')^2 + ...
        norm(R23*R31 - R12.', 'fro')^2 + ...
        norm(R31*R12 - R23.', 'fro')^2;
end