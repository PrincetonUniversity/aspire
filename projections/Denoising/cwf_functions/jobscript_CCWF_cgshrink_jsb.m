function [ denoised_coeff, eigims, egval, num_eigs] = jobscript_CCWF_cgshrink_jsb(index, w_f_rad, CTF_rad_all, basis, sample_points,  mean_coeff, coeff_pos_k, noise_v)

% Main function to estimate covariance (colored noise) and obtain denoised image coefficients by deconvolution (Wiener filter)
% Inputs:
% index: CTF indices for images
% w_f_rad: Whitening filter evaluated along radial quadrature points
% CTF_rad_all: matrix of all radial CTF's, each column is a distinct CTF
% basis: precomputed Fourier Bessel basis
% sample_points: Quadrature points
% mean_coeff: FB coefficients of the mean image
% coeff_pos_k: FB coefficients of images for non-negative angular frequencies
% noise_v: Noise variance
% OUTPUTS:
% denoised_coeff: Cell of FB coeffs of denoised images, each cell corresponds to a non-negative angular frequency
% eigims: Cell of eigenimages corresponding to each angular frequency
% egval: EIgenvalues corresponding to each angular frequency block
% num_eigs: Number of non-zero eigenvalues for each block
% Tejal Jan 2016

nim=length(index);
ang_freqs=basis.ang_freqs;

for i=1:max(index)
    weight(i)=length(find(index==i));
end
denoised_coeff=cell(length(coeff_pos_k), 1);
eigims=cell(length(coeff_pos_k), 1);
egval=cell(length(coeff_pos_k), 1);
num_eigs=0;
for k=unique(ang_freqs)'
    if k>=0
        tmp=coeff_pos_k{k+1};
        clear DD A
        for i=1:max(index)
            DD{i}=real(tmp(:,find(index==i))*tmp(:,find(index==i))'/length(find(index==i)));
            A{i}=calc_fb_CTF(CTF_rad_all(:,i),basis.Phi_ns{k+1}, sample_points);
        end
        
        regu=0;
        if k==0
            disp('Regularizing k=0 block')
            regu=choose_regu_k0(A,weight);
        end
        
        %% Operator shrinkage with regularization for k=0 and KN cutoff
        % [C,relres,iter,num_eig]=solve_cg_shrink_old(A,DD,weight,noise_v, nim,k, regu);
        [C, relres, iter, num_eig]=solve_cg_shrink(A,DD,weight,noise_v, nim,k, regu);
        invW=calc_fb_CTF(w_f_rad.^(-1),basis.Phi_ns{k+1}, sample_points);
        C=invW*C*invW'; toc
        
        %% Ensure positivity of the estimated covariance matrix
        [U,S]=eig(C);
        C=U*diag(max(diag(S),0))*U';
        eigims{k+1}=U;
        egval{k+1}=diag(S);
        num_eigs=num_eigs+num_eig;
        %eigims{k+1}=0; egval{k+1}=0;
        %% Deconvolution
        for i=1:max(index)
            WA{i}=calc_fb_CTF(CTF_rad_all(:,i).*w_f_rad,basis.Phi_ns{k+1}, sample_points);
            h_posk = C*WA{i}'*inv(WA{i}*C*WA{i}'+noise_v*eye(size(WA{i},1)));
            if (k==0)
                denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i)) + repmat(mean_coeff{k+1},1,numel(find(index==i))); %disp('Subtracting mean in deconvolution')
            else
                denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i)) ;
            end
            
        end
    end
end

