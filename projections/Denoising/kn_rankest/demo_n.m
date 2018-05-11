%demo_n

% simple performance curve of rank estimation using the algorithm proposed
% in Kritchman & Nadler, 
% and compared to the MDL-based rank estimator, proposed by Wax & Kailath. 

% YOU CAN CHANGE THE VALUES OF THE FOLLOWING 3 QUANTITIES: 
% p  - number of sensors / dimensionality of each sample
% N  - array of values of number of samples n. 
% lambda_population = % (noise-free) population eigenvalues 

p =16;     
N = [50 100 150 200 300 400 500 600 700 800 1000 1250 1500 1750 2000 2250 2500 3000 3500 4000 4500 5000]; 
lambda_population = [4 2 1 0.25];  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = length(N); 
true_K = length(lambda_population); 

seed = 97289342856; randn('state',seed); % for reproducibility. Change the seed for new results
iter = 5000;  % number of iterations per each value of p,n


sigma = 1;  % noise level
beta = 1;   % beta=1 for real valued observations. Change to beta=2 for complex valued data and noise. 

alpha_KN = 0.1; % confidence level for the algorithm by KN. 

p_WK = []; 
p_KN = []; 

WK  = zeros(L,iter); 
KN = zeros(L,iter); 

for counter = 1:iter

      fprintf('counter %d \n',counter); 
      for i=1:L; 
            n = N(i);
               
            X = zeros(n,p); 
            for j=1:length(lambda_population)  % we put the signals in the first K coordinates
                if beta==1
                    X(:,j)= sqrt(lambda_population(j)) * randn(n,1); 
                else
                    X(:,j)= 1/sqrt(2) * sqrt(lambda_population(j)) * complex(randn(n,1),randn(n,1));  
                end
            end

            if beta==1   % ADD NOISE
                X = X + sigma * randn(n,p); % real noise
            else
                X = X + 1/sqrt(2) * sigma * complex(randn(n,p),randn(n,p)); 
            end

           if p<= n         
                S = 1/n * X' * X;   
                S = (S + S') / 2 ; 
                lambda = eig(S); %eigs(S,'lm',options);            
            else % here p > n better to compute n x n matrix
                S = 1/n * X * X';  
                S = (S + S') / 2 ;             
                lambda = eig(S); %eigs(S,'lm',options);
                lambda(n+1:p) = 0; 
           end
            lambda = sort(lambda,'descend'); 

             KN(i,counter) = KN_rankEst(lambda,n,beta,alpha_KN); 
             WK(i,counter) = rank_est_WK(lambda,p,n); 
 
        
            p_WK(i) = length(find(WK(i,1:counter)==true_K)) / counter; 
            p_KN(i) = length(find(KN(i,1:counter)==true_K)) / counter;            
        end

        figure(12); 
        plot(N,1-p_KN,'rs-',N,1-p_WK,'ko-'); 
        legend('KN','MDL/WK',1);
         xlabel('N','fontsize',16);
         ylabel('Pr(K_{est} \neq  K)','fontsize',16); grid on; 
         title(['\beta= ' num2str(beta) ' K=' num2str(true_K) ' p= ' num2str(p) ', \lambda = [' num2str(lambda_population) ']'],'fontsize',16)
        drawnow; pause(0.001); 
    end

return;     
