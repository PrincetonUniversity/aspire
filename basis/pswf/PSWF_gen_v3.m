function [ PSWF_for_approx, Alpha_Nn, ang_freq, rad_freq ] = PSWF_gen_v3( n_initial, L, beta, phi_approx_err, PSWF_approx_err, x, y )
% This function generates the PSWF Basis functions with positive angular frequencies
% and normalizes them: psi_normalized = c/(2*pi*L)*alpha_{N,n}*psi_{N,n}(x).
% This way the expansion coefficients are given by the standard l2 inner product with
% the equally spaced samples of this basis inside the unit disk.
%   Input: 
%   Output: Psi_Nn: PSWF Basis with order n and s
%           ang_freq: angular frequencies
%           rad_freq: radial frequencies
%
%   Note: The x,y points should be entered as integers between 0 and L instead of 0 and 1(i.e. L times their original values)

r_2d_grid_on_the_circle = sqrt(x.^2+y.^2)/L;
theta_2d_grid_on_the_circle =  angle(x/L + 1i*y/L);

alpha_all = [];
PSWF_for_approx = [];

c = beta * pi*L;

N = 0;
n = n_initial;
n_order_length_vec = [];
while(1)
    [PSWF , alpha] = PSWF_2D(N,n,c,r_2d_grid_on_the_circle,theta_2d_grid_on_the_circle,phi_approx_err);
    
    lambda = (c/(2*pi))^2 * abs(alpha).^2;
    gamma = sqrt(abs(lambda./(1-lambda))); 

    n_end = find(gamma<=PSWF_approx_err,1,'first');
    
    if (n_end < 2)
        break;
    end
    
    if (~isempty(n_end))    
        n_order_length_vec = [n_order_length_vec,n_end-1];
        alpha_all = [alpha_all;alpha];
        PSWF_for_approx = [PSWF_for_approx,PSWF];        
        N = N +1;
        % clc; 
        % display(['Generating PSWFs for angular index: ',num2str(N)]);
        n = n_end;
    else
        % The case where the initial n isn't large enough -> Try again.        
        n = 2*n;
    end
end

lambda = (c/(2*pi))^2 * abs(alpha_all).^2;
gamma = sqrt(abs(lambda./(1-lambda))); 
    
PSWF_for_approx = PSWF_for_approx(:,gamma>PSWF_approx_err);
alpha_for_approx = alpha_all(gamma>PSWF_approx_err);

PSWF_for_approx = bsxfun(@times,PSWF_for_approx,(beta/2)*alpha_for_approx.'); 
Alpha_Nn = alpha_for_approx;

ang_freq = [];
rad_freq = [];
N = 0;
for i = 1:length(n_order_length_vec)
    ang_freq = [ang_freq; N*ones(n_order_length_vec(i),1)];
    rad_freq = [rad_freq; (1:n_order_length_vec(i)).'];
    N = N+1;
end
 
