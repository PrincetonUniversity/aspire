function [phiZeros] = findInitialNodes(x,n,c,phi_approx_err,idxForQuadNodes)
%% Generate points for numerical integration
N = 0;

%% Calculate approx_len
% - Definition for the approximation function of the d's decay
d_decay_approx_fun = @(N,n,c,j) c^2./(16*(j.^2 + j*(2*n+N+1)) - c^2);

% - Finding the first point of the decay according to the approximation
first_idx_for_decrease = ceil((((2*n+N+1)^2+c^2/2)^0.5 - (2*n+N+1))/2);     % Analytic solution for d_decay_approx_fun<1 for every j>first_idx_for_decrease

% - Finding the number of points needed to achieve the chosen error
d_approx = d_decay_approx_fun(N,n,c,first_idx_for_decrease);
d_decay_idx_counter = first_idx_for_decrease;

while(d_approx>phi_approx_err)
    d_decay_idx_counter = d_decay_idx_counter + 1;
    d_approx = d_approx * d_decay_approx_fun(N,n,c,d_decay_idx_counter);
end

approx_len = n + 1 + d_decay_idx_counter;

%% Funcion definitions
Pn = @(n,alpha,beta,x) j_polynomial(length(x),n,alpha,beta,x);

%% Generate Matrix for calculation of d_k coefficients
B_N = Generate_BN_mat(N,c,approx_len);
[d_vec, xi] = eig(B_N);

%% Sort eigenvalues by descending order
[val,idx] = sort(diag(xi),'descend');
d_vec = d_vec(:,idx);
for i=1:approx_len
    xi(i,i) = val(i);
end

%% Generate T_n basis functions on grid points ,and the derivatives of T_n.
j=0:(approx_len-1);

% - Calculation of functions
T_x_mat = @(x,N,j,approx_len) ((x.^(N+1/2))*((2*(2*j+N+1)).^(1/2))) .* Pn(approx_len-1,N,0,1-2*x.^2);

%% Find zero's of phi with zeros equal to number of quadrature nodes (exluding x=0)
phiForQuadNodes = @(x) ( T_x_mat(x,N,j,approx_len) * d_vec(:,idxForQuadNodes) );

% - Find K+1 distinct zeros
funVec = phiForQuadNodes(x);
signFlippingVec = find(sign(funVec(1:end-1)) ~= sign(funVec(2:end))).';
phiZeros = [];
for i = signFlippingVec
%     i/numel(x)
    newZero = x(i)-phiForQuadNodes(x(i))*( (x(i+1)-x(i))/(phiForQuadNodes(x(i+1))-phiForQuadNodes(x(i))) );
    phiZeros = [phiZeros newZero];
    if (numel(phiZeros)==(idxForQuadNodes-1))
        break;
    end
end
phiZeros = phiZeros(1:idxForQuadNodes-1);

% figure; plot(x,phiForQuadNodes(x)); grid on
% hold on; scatter(phiZeros,zeros(1,numQuadPts),[],[1 0 0],'x')
end