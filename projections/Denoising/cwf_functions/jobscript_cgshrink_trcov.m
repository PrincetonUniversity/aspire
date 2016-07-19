function [ denoised_coeff, err_cov, mse_tot, cov, err_rhs, iters, res, rhsnt ] = jobscript_cgshrink_trcov(index, CTF_rad_all, basis, sample_points,  mean_coeff, coeff_pos_k, tr_cov, noise_v)

% Main function to estimate covariance (white noise) and obtain denoised image coefficients by deconvolution (Wiener filter)
% Also calculate error in estimated covariance
%
% Inputs:
% index: CTF indices for images
% w_f_rad: Whitening filter evaluated along radial quadrature points
% CTF_rad_all: matrix of all radial CTF's, each column is a distinct CTF
% basis: precomputed Fourier Bessel basis
% sample_points: Quadrature points
% mean_coeff: FB coefficients of the mean image
% coeff_pos_k: FB coefficients of images for non-negative angular frequencies
% noise_v: Noise variance
%
% OUTPUTS:
% denoised_coeff: Cell of FB coeffs of denoised images, each cell corresponds to a non-negative angular frequency
% err_cov: (Relative) Error in estimated covariance (block-wise)
% mse_tot: (Relative) Error in estimated covariance - total. The measure is the
% squared Frobenius norm
% cov: Estimated covariance (blocks)
% err_rhs: Error in the RHS of the linear system solved. See eqn. 9 in http://arxiv.org/pdf/1602.06632v3.pdf
% iters: Number of CG iterations
% res: Residual in CG
% rhsnt: Normalized RHS
% Tejal Jan 2016

nim=length(index);
ang_freqs=basis.ang_freqs;

for i=1:max(index)
    weight(i)=length(find(index==i));
end

denoised_coeff=cell(length(coeff_pos_k), 1);
cov=cell(length(coeff_pos_k), 1);
err_tot=0;
ctr_tot=0;
eigims=cell(length(coeff_pos_k), 1);

for k=unique(ang_freqs)'
    
    if k>=0
        tmp=coeff_pos_k{k+1};
        tr_cov_k=tr_cov{k+1};
        clear DD A
        for i=1:max(index)
            
            DD{i}=real(tmp(:,find(index==i))*tmp(:,find(index==i))'/length(find(index==i)));
            [UDD,SDD]=eig((DD{i}));
            
            % Shrinkage for sample covariance
            n_i=length(find(index==i));
            if (k==0)
                gamma=size(DD{i},1)/n_i;
            else
                gamma=size(DD{i},1)/(2*n_i);
            end
            
            A{i}=calc_fb_CTF(CTF_rad_all(:,i),basis.Phi_ns{k+1}, sample_points);
        end
        
        regu=0;
        %            if k==0
        %               disp('Regularizing k=0 block')
        %               regu=choose_regu_k0(A,weight)
        %            end
        %No regularization for white noise case
        %regu=0;
        
        [C, e_rhs, iter, relres, rhsn]=solve_cg_shrink_deb(A,DD,tr_cov{k+1},weight,noise_v, nim,k,regu);  % Changed, solve_LS_c from solve_LS which doesn't include weights
        
        
        %% Ensure positivity of the estimated covariance matrix
        [U,S]=eig(C);
        C=U*diag(max(diag(S),0))*U';
        cov{k+1}=C;
        err_rhs(k+1)=e_rhs; % Normalized error should show 1/sqrt(n) by LLN
        err_cov(k+1)=norm(C-tr_cov_k,'fro').^2/norm(tr_cov_k,'fro').^2;
        err_tot=err_tot+norm(C-tr_cov_k,'fro').^2;
        ctr_tot=ctr_tot+norm(tr_cov_k,'fro').^2;
        iters(k+1)=iter;
        res{k+1}=relres;
        rhsnt(k+1)=rhsn;
        for i=1:max(index)
            
            h_posk = C*A{i}'*inv(A{i}*C*A{i}'+noise_v*eye(size(A{i},1)));
            
            if (k==0)
                denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i)) + repmat(mean_coeff{k+1},1,numel(find(index==i))); %disp('Subtracting mean in deconvolution')
            else
                denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i)) ; %disp('Subtracting mean in deconvolution')
            end
            
        end
    end
end

mse_tot=err_tot/ctr_tot;
