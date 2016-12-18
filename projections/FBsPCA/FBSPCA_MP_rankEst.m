function [ UU, Freqs, Rad_Freqs, W ] = FBSPCA_MP_rankEst( P, U, D, freqs, rad_freqs, nv)
%This function select signal components
%   Input: P: number of images
%          U: principal components 
%          D: eigenvalues
%          freqs: associated angular frequencies
%          rad_freqs: associated radial frequencies
%          nv: noise variance sigma^2.
%   Output: UU: Eigen functions
%           Freqs: associated angular frequency
%           Rad_Freqs: associated radial frquency 
%           W: Linear filter weights
%Zhizhen Zhao Sep 2012

%%%For zero angular frequency
K = MP_rankEst(D(freqs==0), P, nv); %estimate number of components to keep
gamma=length(D(freqs==0))/P;
freqs_tmp=freqs((freqs==0), :);
rad_freqs_tmp=rad_freqs((freqs==0));
D_tmp=D(freqs==0);
U_tmp=U(:, freqs==0);
freqs_tmp=freqs_tmp(1:K);
rad_freqs_tmp=rad_freqs_tmp(1:K);
D_tmp=D_tmp(1:K);
% Estimate the lambda_k
l_k=0.5*((D_tmp-(gamma+1)*nv)+sqrt(((gamma+1)*nv-D_tmp).^2-4*gamma*nv^2));
% SNR_k
SNR_k=l_k/nv; 
% SNR_{k, \gamma}
SNR=(SNR_k.^2-gamma)./(SNR_k+gamma); 
U_tmp=U_tmp(:, 1:K);
weight=1./(1+1./SNR);
Freqs=freqs_tmp;
Rad_Freqs=rad_freqs_tmp;
UU=U_tmp;
W=weight;
%For non-zero angular frequencies
for i=1:max(freqs)
    [ K ]=MP_rankEst(D(freqs==i), 2*P, nv);
    gamma=length(D(freqs==i))/(2*P);
    if K~=0    
        freqs_tmp=freqs((freqs==i));
        rad_freqs_tmp=rad_freqs((freqs==i));
        D_tmp=D(freqs==i);
        U_tmp=U(:, freqs==i );
        freqs_tmp=freqs_tmp(1:K);
        rad_freqs_tmp=rad_freqs_tmp(1:K);
        D_tmp=D_tmp(1:K);
        % Estimate the lambda_k
        l_k=0.5*((D_tmp-(gamma+1)*nv)+sqrt(((gamma+1)*nv-D_tmp).^2-4*gamma*nv^2));
        % SNR_k
        SNR_k=l_k/nv;
        % SNR_{k, gamma}
        SNR=(SNR_k.^2-gamma)./(SNR_k+gamma);
        U_tmp=U_tmp(:, 1:K);
        weight=1./(1+1./SNR);
    else
        continue;
    end;
    Freqs = [Freqs; freqs_tmp];
    Rad_Freqs = [Rad_Freqs; rad_freqs_tmp];
    UU = [UU, U_tmp];
    W=[W; weight];
end;
      
end

