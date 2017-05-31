function [ denoised_coeff, skip_flags]  = jobscript_CCWF_cgshrink_test(index, CTF_rad_all, basis, sample_points,  mean_coeff,  noise_v, nbatch, nim, fpath)
%Tejal Bhamre, Oct 2015
%INPUTS:
%index : CTF indices
%CTF_rad_all: CTF evaluated along radial quadrature points
%basis: Precomputed basis functions
%sample_points: Quadrature points
%mean_coeff: Coefficients of mean image
%noise_v: Noise variance
%fpath: Destination of folder for output files
%OUTPUTS:
%denoised_coeff: Coeffs of denoised images
%skip_flags: Flags for skipping conjugate gradient


ang_freqs=basis.ang_freqs;

for i=1:max(index)
    weight(i)=length(find(index==i));
end
%disp('Change nim for true cov here later')
%denoised_coeff=cell(length(coeff_pos_k), 1);
%eigims=cell(length(coeff_pos_k), 1);
% Only pos k, to be compatibel with Jane's ffb
for k=unique(ang_freqs)'
    
    if k>=0
        filename=fullfile(fpath, sprintf('coeff_pos%d',k));
 	load(filename);
        tmp=coeff_pos_k_new;
        %tr_c=tr_coeff_pos_k{k+1};
        %tr_cov_k=real(tr_c*tr_c')/(1e5);
       % tr_sample_k=tr_samplecov{k+1};
        %tmp1=conj(tmp);
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


           % tr_DD{i}=real(tr_sample_k(:,find(index==i))*tr_sample_k(:,find(index==i))'/length(find(index==i)));
           
            A{i}=calc_fb_CTF(CTF_rad_all(:,i),basis.Phi_ns{k+1}, sample_points);
        end
        
            regu=0; 
%            if k==0
%               disp('Regularizing k=0 block')
%               regu=choose_regu_k0(A,weight)
%            end
 %No regularization for white noise case
            regu=0;

	    [C, relres, iter, num_eig, skip_flag]=solve_cg_shrink(A,DD,weight,noise_v, nim,k,regu);  % Changed, solve_LS_c from solve_LS which doesn't include weights

   
         %% Ensure positivity of the estimated covariance matrix
        [U,S]=eig(C);
        C=U*diag(max(diag(S),0))*U';
	skip_flags(k+1)=skip_flag;
	%cov{k+1}=C;
        %eigims{k+1}=U;        
        %err_cov(k+1)=norm(C-tr_cov_k,'fro').^2/norm(tr_cov_k,'fro').^2;
        %err_tot=err_tot+norm(C-tr_cov_k,'fro').^2;
        %ctr_tot=ctr_tot+norm(tr_cov_k,'fro').^2;

        for i=1:max(index)

            h_posk = C*A{i}'*inv(A{i}*C*A{i}'+noise_v*eye(size(A{i},1)));
            
            if (k==0)
            denoised_coeff(:,find(index==i)) = h_posk*tmp(:,find(index==i)) + repmat(mean_coeff{k+1},1,numel(find(index==i))); %disp('Subtracting mean in deconvolution')
            else
            denoised_coeff(:,find(index==i)) = h_posk*tmp(:,find(index==i)) ; %disp('Subtracting mean in deconvolution')
            end
            % Fix this
        end
	    filename=fullfile(fpath, sprintf('den_coeff_pos%d',k));
            save(filename,'denoised_coeff','-v7.3');
	    clear denoised_coeff
    end
end

%mse_tot=err_tot/ctr_tot;
