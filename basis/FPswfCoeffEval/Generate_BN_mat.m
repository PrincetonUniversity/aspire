function [ B_N ] = Generate_BN_mat(N,c,approx_len)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Defenitions %%
% my_factorial = @(n,k) factorial(n)./(factorial(n-k).*factorial(k));
% h = @(N,n) (2*(2*n+N+1).*(nchoosek(n+N,n)).^2).^(1/2);
% h = @(N,n) (2*(2*n+N+1).*(my_factorial(n+N,n)).^2).^(1/2);
h = @(N,n) (2*(2*n+N+1)).^(1/2);

k = @(N,n) (N+2*n+1/2).*(N+2*n+3/2);

% gamma_plus_1 = @(N,n) -(((n+N+1).^2).*h(N,n))./((2*n+N+1).*(2*n+N+2).*h(N,n+1));
gamma_plus_1 = @(N,n) -((((n+N+1).^2).*h(N,n))./((2*n+N+1).*(2*n+N+2).*h(N,n+1))).*((n+1)./(n+N+1));    % Simplified version
% gamma_plus_1 = @(N,n) -((((n+N+1).^2))./((2*n+N+1).*(2*n+N+2)));    % Original Slepian paper version
if (N==0)
    gamma_0 = @(N,n) 1/2;
else
    gamma_0 = @(N,n) (2*n.*(n+1) + N*(2*n+N+1))./((2*n+N).*(2*n+N+2));
end
% gamma_minus_1 = @(N,n) -((n.^2).*h(N,n))./((2*n+N+1).*(2*n+N).*h(N,n-1));     % Simplified version
gamma_minus_1 = @(N,n) -(((n.^2).*h(N,n))./((2*n+N+1).*(2*n+N).*h(N,n-1))).*((n+N)./n);
% gamma_minus_1 = @(N,n) -(((n.^2))./((2*n+N+1).*(2*n+N)));           % Original Slepian paper version

b_N_above_diag = @(N,n,c) -c^2*gamma_minus_1(N,n);
b_N_on_diag = @(N,n,c) -(k(N,n) + c^2*gamma_0(N,n));
b_N_below_diag = @(N,n,c) -c^2*gamma_plus_1(N,n);

%% Generate B_N matrix %%

B_N = zeros(approx_len);

B_N(1,1) = b_N_on_diag(N,0,c);

i = 2:approx_len;
b_N_above_diag_vec = b_N_above_diag(N,i-1,c);
b_N_on_diag_vec = b_N_on_diag(N,i-1,c);
b_N_below_diag_vec = b_N_below_diag(N,i-2,c);

for i=2:approx_len
    B_N(i-1,i) = b_N_above_diag_vec(i-1);
    B_N(i,i) = b_N_on_diag_vec(i-1);
    B_N(i,i-1) = b_N_below_diag_vec(i-1);
end

end

