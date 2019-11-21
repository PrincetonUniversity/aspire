function [ func_2d, T] = SOFT_yuan( A, B, L, bandl)
% Fourier transform coefficients to function values over SO(3)
% Input:
%   L: maximum bandlimit, starting from 0 to L-1
%   coeff: vector of coefficients to be transformed

% Output:
%   func: function values for uniformly spaced euler angles, where 
%   alpha_i=i*pi/L, beta_k=pi*(2*k+1)/(4*L), gamma_j=j*pi/L, and 
%   i,j,k=[0,...,2*L-1], and j iterates the fastest, i, the second, and k
%   the slowest.
format long;
N = 2*bandl-1;

% direct summation for test O(L^6)
% func = zeros(1,(2*L)^3);
% iter = 0;
% 
% %                 % why this is not working?
% %                 b = conj(B(:,(l1-1)^2+1:l1^2))';
% %                 a = A(:,(l1-1)^2+1:l1^2);
% %                 fmmp = b*a;
%                 %wigd = wignerd(l1-1,-beta_k);
% %                 m = -(l1-1):l1-1;
% %                 func(iter) = func(iter)+exp(-1i*m*gamma_j2)*(wigd.*fmmp)*...
% %                    exp(-1i*m*alpha_j1)'; % be careful of m and mp indices
% tic
% for k = 0:N
%     beta_k = pi*(2*k+1)/(4*L);
%     for j1 = 0:N
%         alpha_j1 = pi*j1/L;
%         for j2 = 0:N
%             gamma_j2 = pi*j2/L;
%             iter = iter+1;
%             for l1 = 1:L 
%                 [~,wigd,~] = wignerD(l1-1,0,beta_k,0);
%                 for m = -(l1-1):l1-1 % good with this naive way
%                     for mp = -(l1-1):l1-1
%                         b = conj(B(:,(l1-1)^2+m+l1))';
%                         a = A(:,(l1-1)^2+mp+l1);
%                         fmmp = b*a;
%                         func(iter) = func(iter)+exp(-1i*m*gamma_j2)*wigd(m+l1,mp+l1)*...
%                             exp(-1i*mp*alpha_j1)*fmmp;
%                     end
%                 end
%             end
%         end
%     end
% end
% t_direct = toc
% func = conj(func);
%max(abs(func-func_sb)) %10^-21 checked
 
% a little faster according to Roy with 2d fft 
f_beta = zeros(2*L-1,2*L-1,2*bandl); % row index m, column index mp
func_2d = zeros(2*bandl*2*bandl,2*bandl);
alpha = pi*(0:N)/bandl;
[gamma,alpha] = ndgrid(alpha,alpha);
tic
for k = 0:N
    beta_k = pi*(2*k+1)/(4*bandl);
    for l1 = 1:L
        %[~,wigd,~] = wignerD(l1-1,0,beta_k,0); % use easy spin
        wigd = wignerd(l1-1,beta_k);
        for m = -(l1-1):l1-1 % good with this naive way
            for mp = -(l1-1):l1-1
                b = conj(B(:,(l1-1)^2+m+l1))';
                a = A(:,(l1-1)^2+mp+l1);
                fmmp = b*a;
                f_beta(m+L,mp+L,k+1) = f_beta(m+L,mp+L,k+1)+fmmp*wigd(m+l1,mp+l1);
            end
        end
    end
    func_2d(:,k+1) = nufft2d2((2*bandl)^2,gamma(:),alpha(:),-1,1e-15,2*L-1,2*L-1,f_beta(:,:,k+1));
end    
t_soft = toc
func_2d = func_2d(:); % good to 10^21 now
%err_soft = max(abs(func_2d'-func))

% % matrix multiplication
% func_sb = zeros(1,(2*L)^3);
% iter = 0;
% tic
% for k = 0:N
%     beta_k = pi*(2*k+1)/(4*L);
%     for j1 = 0:N
%         alpha_j1 = pi*j1/L;
%         for j2 = 0:N
%             gamma_j2 = pi*j2/L;
%             iter = iter+1;
%             R = euler2rotationMatrix(alpha_j1, beta_k, gamma_j2, 'zyz');
%             func_sb(iter) = trace(conj(B(:,1:L^2))*(getSHrotMtx(R,L-1,'complex'))*(A(:,1:L^2))');
%         end
%     end
% end
% t_matrix = toc
% err_matrix = max(abs(func_sb-func_2d'))

T = t_soft;
%E = err_matrix;

end

% k = 0;
% j1=1;
% j2 =1;
% beta_k = pi*(2*k+1)/(4*L);
% alpha_j1 = pi*j1/L;
% gamma_j2 = pi*j2/L;
% func_i = nufft2d2(1,gamma_j2,alpha_j1,-1,1e-15,2*L-1,2*L-1,f_beta(:,:,k+1)); % why minus sign?
% 
% k1=-L+1:L-1;
% k2=k1;
% [K1,K2]=ndgrid(k1,k2);
% k1=K1(:); clear K1;
% k2=K2(:); clear K2;
% cg3=zeros((2*L)^2,1);
% gamma = gamma(:);
% alpha = alpha(:);
% tic
% for j=1:(2*L)^2
% 	x1j=gamma(j);
% 	x2j=alpha(j);
% 	cg3(j)=sum(reshape(f_beta(:,:,k+1),(2*L-1)^2,1).*exp(-1i*(k1*x1j+k2*x2j)) );
% end 
% toc
% figure()
% subplot(2,1,1)
% plot(1:64,cg3);
% subplot(2,1,2)
% plot(1:64,func_2d(1:64));

% more faster with SOFT 1d fft O(L^4)

%end

