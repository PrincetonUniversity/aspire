function [ PSWF_for_approx, Alpha_Nn, ang_freq, rad_freq ] = PSWF_2D_gen_with_truncation( L, beta, phi_approx_err, truncation_param, x, y )
%This function generates the PSWF basis with positive angular frequencies
%for those indices which correspont to the truncation rule.
%   Input: 
%   Output: Phi_ns: PSWF Basis with order n and s
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
n = L+1;
n_order_length_vec = [];
while(1)
    [PSWF , alpha] = PSWF_2D(N,n,c,r_2d_grid_on_the_circle,theta_2d_grid_on_the_circle,phi_approx_err);
    
    lambda = (c/(2*pi))^2 * abs(alpha).^2;
    gamma = sqrt(abs(lambda./(1-lambda))); 

    n_end = find(gamma<=truncation_param,1,'first');
    
    if (n_end < 2)
        break;
    end
    
    if (~isempty(n_end))    
        n_order_length_vec = [n_order_length_vec,n_end-1];
        alpha_all = [alpha_all;alpha];
        PSWF_for_approx = [PSWF_for_approx,PSWF];        
        N = N +1
        n = n_end;
    else
        % The case where the initial n isn't large enough -> Try again.        
        n = 2*n;
    end
end

lambda = (c/(2*pi))^2 * abs(alpha_all).^2;
gamma = sqrt(abs(lambda./(1-lambda))); 
    
PSWF_for_approx = PSWF_for_approx(:,gamma>truncation_param);
alpha_for_approx = alpha_all(gamma>truncation_param);

Alpha_Nn = alpha_for_approx;

ang_freq = [];
rad_freq = [];
N = 0;
for i = 1:length(n_order_length_vec)
    ang_freq = [ang_freq; N*ones(n_order_length_vec(i),1)];
    rad_freq = [rad_freq; (1:n_order_length_vec(i)).'];
    N = N+1;
end
 