function results=test_cryo_clmatrix_ml_gpu(shiftslist,KNNs,matname)
%
% TEST_CRYO_CLMATRIX_ML     Test performance of cryo_clmatrix_ml
%
% Note that the function tests the GPU implementaion.
%
% results=test_cryo_clmatrix_ml(shiftslist,KNNs,matname)
%   Compare the performance of maximum likelihood common lines detection to
%   correlation-based detection. The maximimum likelihood uses KNNs
%   neighbors for updating the ML estimates. The tests images are shifts
%   according to the shifts in the array shiftslist.
%   Since the running time of this function is rather long, intermediate
%   ourpur is saved to the MAT file matname. Default name is
%   test_cryo_clmatrix_ml.mat
%
%   Default shiftslist is [0,5]. Pass an empty array [] to use default
%   shifts. Default KNNs is 200.
%
%   Returns MATLAB table result, which can be displayed using 
%           struct2table(results)
%
%   The function generates 300 images of size 129x129, add noise of varying
%   SNRand compares the detection rate of correlation based detection and
%   of maximum likelihood-based detection. The noisy projections are masked
%   before searching for common lines.
%   
%   Examples:
%   results=test_cryo_clmatrix_ml;          % shifts are [0,5], KNNs=200
%   results=test_cryo_clmatrix_ml([0 10]);  % shifts are [0,10], KNNs=200
%   results=test_cryo_clmatrix_ml([],100);  % shifts are [0,5], KNNs=100
%
% Yoel Shkolnisky, July 2015.


if (nargin==0) || (nargin>0 && isempty(shiftslist))
    shiftslist=[0,5];
end

if nargin<2
    KNNs=[200]; % In the future, just put multiple values in this list 
                   % to test several KNN values.
end

if nargin<3
    matname='test_cryo_clmatrix_ml_gpu';
end

% Fixed parameters of the experiments
n=129;  % Size of each image
K=300;  % Number of images
shift_step_2d=1;    % Shift step resolution to shift the input images 
                    % (in pixels).

SNRlist=[1/4,1/8,1/16,1/32]; % SNR values to use.
testidx=1;
Ntests=numel(SNRlist)*numel(shiftslist)*numel(KNNs); % Total number of test that will be executed.

results=struct; % Fields are snr,maxshift,KNN,rate_corr,rate_ml,t_corr,t_ml
for snridx=1:numel(SNRlist)
    for maxshiftidx=1:numel(shiftslist);
        for KNNidx=1:numel(KNNs);
            
                SNR=SNRlist(snridx);
                max_shift_2d=shiftslist(maxshiftidx);
                M=ceil(1/(2*(sind(1/4))^2)); % Number of underying clean signals used.

                log_message('*********************************');               
                log_message('Starting test %d/%d',testidx,Ntests);
                log_message('\t SNR=%d',SNR);
                log_message('\t max_shift_2d=%d',max_shift_2d);


                %% Generate projections
                initstate; % So we get the same results every time for reproducibility.
                rots_ref = rand_rots(K);  % Generate random uniform rotations.
                
                load cleanrib
                volref=real(volref);
                volref=Downsample(volref,[n,n,n]);
                projs=cryo_project(volref,rots_ref,n);
                projs=permute(projs,[2 1 3]); % Swap dimensions for compitability with old gen_projections.
                [projs,~]=cryo_addshifts(projs,[],max_shift_2d,shift_step_2d);
                noisy_projs=cryo_addnoise(projs,SNR,'gaussian');
                masked_projs=mask_fuzzy(noisy_projs,floor(n/2)); % Applly circular mask
                
                % Compute polar Fourier transform, using radial resolution n_r and angular
                % resolution n_theta. n_theta is the same as above.
                n_theta=360;
                n_r=round(n/2);
                [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
                
                % Find common lines from projections
                max_shift_1d=ceil(2*sqrt(2)*max_shift_2d);
                shift_step_1d=shift_step_2d;

                tic;
                clstack = cryo_clmatrix_gpu(pf,K,1,max_shift_1d,shift_step_1d);
                t_corr=toc;
                
                % Find reference common lines and compare
                [ref_clstack,~]=clmatrix_cheat(rots_ref,n_theta);
                prop_corr=comparecl( clstack, ref_clstack, n_theta, 5 );
                log_message('Percentage of correct common lines: %f%%',prop_corr*100);

                KNN=KNNs(KNNidx);
                max_itr=6;
                verbose=1;
                
                [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections                
                psd=cryo_noise_estimation_pfs(noisy_projs,n_r,n_theta);
                tic;
                [clstack2]=...
                    cryo_clmatrix_ml_gpu(pf,M,KNN,diag(psd),max_itr,verbose,max_shift_1d,shift_step_1d);                
                t_ml=toc;
                prop_ml=comparecl( clstack2, ref_clstack, n_theta, 5 );
                log_message('Percentage of correct common lines: %f%%',prop_ml*100);
                
                results(testidx).snr=SNR;
                results(testidx).maxshift=max_shift_2d;
                results(testidx).KNN=KNN;
                results(testidx).rate_corr=prop_corr;
                results(testidx).rate_ml=prop_ml;
                results(testidx).t_corr=t_corr;                
                results(testidx).t_ml=t_ml;

                log_message('\t KNN=%d',KNN);
                log_message('\t rate_corr=%4.2f',prop_corr*100);
                log_message('\t rate_ml=%4.2f',prop_ml*100);
                log_message('\t t_corr=%d',t_corr);
                log_message('\t t_ml=%d',t_ml);
                
                save(matname,'results');
                
                testidx=testidx+1;
                
                log_flush;

        end
    end
end
