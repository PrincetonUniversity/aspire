function [ PSWF_for_approx, Alpha_Nn, ang_freq, rad_freq ] = PSWF_2D_radial_quad( L, beta, phi_approx_err, truncation_param, radialQuadPts )
%This function generates the normalized radial part of the PSWFs for the fast coefficient ealuation method.
%   Input: 
%   Output: Phi_ns: PSWF Basis with order n and s
%           ang_freq: angular frequencies
%           rad_freq: radial frequencies

alpha_all = [];
PSWF_for_approx = [];

c = beta * pi*L;

N = 0;
n = L+1;
n_order_length_vec = [];
while(1)
    [PSWF , alpha] = PSWF_2D(N,n,c,radialQuadPts,zeros(numel(radialQuadPts),1),phi_approx_err);
    
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
        N = N +1;
%         clc; 
%         display(num2str(N));
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
 