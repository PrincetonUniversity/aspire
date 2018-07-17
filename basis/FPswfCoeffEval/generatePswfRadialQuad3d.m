function [quadRulePts,quadRuleWts] = generatePswfRadialQuad3d(n,c,phi_approx_err,lambdaMax)
%% Generate points for numerical integration
[x,w]=lgwt(20*n,0,1);    % Legendre-Gauss nodes for numerical integration   
x = flipud(x);
w = flipud(w);

N = 0.5;

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

% - Calculation of derivatives
T_x_derivative_mat = ((x.^(N+3/2))*(-2*(N+j+1).*((2*(2*j+N+1)).^(1/2)))).*[zeros(length(x),1),Pn(approx_len-2,N+1,1,1-2*x.^2)] + (N+1/2)*((x.^(N-1/2))*((2*(2*j+N+1)).^(1/2))).*Pn(approx_len-1,N,0,1-2*x.^2);

%% Calculate phi - Generalized PSWF's ,and the derivatives of phi.
phi = T_x_mat(x,N,j,approx_len) * d_vec(:,1:(n+1));
phi_derivatives = T_x_derivative_mat * d_vec(:,1:(n+1));

%% Calculate original problem eigenvalues gamma_N,n
% - First (largest) eigenvalue numerical calculation
[~,max_phi_x_idx] = max(abs(phi(:,1)));
x_for_calc = x(max_phi_x_idx);

K_operator = @(nu,x) besselj(nu,x).*sqrt(x);
Right_hand_side_integral = (w.*K_operator(N,c*x_for_calc*x)).' * phi(:,1);
lambda_N_1 = Right_hand_side_integral/phi(max_phi_x_idx,1);

% - Calculate the rest of the eigenvalues based on the first
upper_integral_values = diag( (((w.*x)*ones(1,n)).*phi_derivatives(:,1:(end-1))).' * phi(:,2:end) );
lower_integral_values = diag( (((w.*x)*ones(1,n)).*phi(:,1:(end-1))).' * phi_derivatives(:,2:end) );
lambda_N = [lambda_N_1 ; lambda_N_1 * cumprod(upper_integral_values./lower_integral_values)];

alpha_N = lambda_N*2*pi*((1i)^N) / sqrt(c);
% alpha_N = lambda_N*2*pi*((1i)^N) / c;

%% Test orthonormality of phi's
% Inner_prod_matrix = ((w*ones(1,n+1)).*phi).' * phi;
% figure; bar(abs(eig(Inner_prod_matrix)));

%% Find indices of suffient radial functions
% K = find(sqrt(c)/2/pi * abs(alpha_N) < lambdaMax,1,'first');
K = find(c/2/pi * abs(alpha_N) < lambdaMax,1,'first');
if (mod(K,2)==0)
    K=K+1;
end

idxForQuadNodes = (K+1)/2;
numQuadPts = idxForQuadNodes - 1;

%% Find inital quadrature nodes - zeros of phi^(c/2)_((K+1)/2)
phiZeros = findInitialNodes(x,n,c/2,phi_approx_err,idxForQuadNodes);

%% Solve system Ax=b (in LS sense) for initial quadrature weights
phiForQuadweights = @(x) ( T_x_mat(x,N,j,approx_len) * d_vec(:,1:K-1) );
% - Compute integrals on right hand-side of equation numerically
b = phiForQuadweights(x).' * (w.*x);

% - Compute matrix A from left hand-side of equation using initial quadrature nodes
A = bsxfun(@times,phiForQuadweights(phiZeros.').',(phiZeros));
initQuadWeights = A\b;

tol = 1e-16;
objFun = @(quadRule) ( bsxfun(@times,phiForQuadweights(quadRule(1:numQuadPts)).',(quadRule(1:numQuadPts)).') * quadRule(numQuadPts+1:end) - b);
% options = optimoptions('lsqnonlin','algorithm','levenberg-marquardt','display','iter-detailed','MaxIter',1e3,'MaxFunEvals',inf,'Tolx',1e-16,'TolFun',1e-16);
options = optimoptions('lsqnonlin','display','iter-detailed','MaxIter',1e3,'MaxFunEvals',inf,'TolFun',tol,'TolX',tol);
quadRuleFinal = lsqnonlin(objFun,[phiZeros.'; initQuadWeights],[],[],options);
quadRulePts = quadRuleFinal(1:numQuadPts);
quadRuleWts = quadRuleFinal(numQuadPts+1:end);
end