function g=nufft_t_3d(beta,x,precision)
%
% Compute the adjoint non-equally spaced FFT in three dimensions.
%
% The function computes the sums
%               n/2
%        g(j) = sum beta(k1,k2,k3)*exp(i*(k1,k2,k3)*x(j))
%             k1,k2,k3=-n/2
% for j=1,...,n.
%
% The complexity of the algorithm is O(n^3*log(1/eps)+(m*n)^3(log(n)))
% The code is implemented such that eps corresponds (by default) to single
% precision. m is the oversampling factor (m=2).
%
% Input parameters:
%    beta       Volume to transform - real or complex. beta is assumed to
%               be a cubed volume.
%    x          Sampling points. Real array of size Nx3 of numbers in the
%               range [-pi,pi]^3. 
%    precision  'double' or 'single'. Default is 'single'.
%    
%    beta and x need not be the same length, but x should be longer than
%    beta.
%
% Output:
%    f        The sums defined above.
%
% Yoel Shkolnisky, January 2009.
%
% Revisions:
% Filename changed from nufft_t3_v2 to nufft_t_3d. (Y.S. December 22, 2009)
%

if nargin<3
    precision='single';
end

if (~strcmpi(precision,'single')) && (~strcmpi(precision,'double'))
    precision='single';
    warning('MATLAB:unkownOption','Unrecognized precsion. Using ''single''.');
end

if strcmpi(precision,'double')
    b=1.5629;
    m=2;
    q=28;
else %single precision
    b=0.5993;
    m=2;
    q=10;
end


if mod(q,2)==1
    error('Set q to an even integer');
end

if ndims(beta)~=3
    error('beta must be a cubed volume');
end

if (size(beta,1)~=size(beta,2)) || (size(beta,1)~=size(beta,2))
    error('beta must be a cubed volume');
end

if length(size(x))~=2
    error('x must be a mx3 array')
end

if size(x,2)~=3
    error('x must be a mx3 array');
end

n=size(beta,1);
len_x=size(x,1);

low_idx=-ceil((n-1)/2);
high_idx=floor((n-1)/2);
idx=low_idx:high_idx;

% Q=zeros(len_x,q+1,q+1,q+1);
g=zeros(len_x,1);
nu=round(x*m*n/(2*pi));

% Do not precompute Q -- it is too big.
% c1=(2*sqrt(b*pi))^3;
% xx=x.*m.*n/(2*pi)-nu;
% kk=zeros(len_x,3);
% for k1=-q/2:q/2
%     kk(:,1)=k1;
%     for k2=-q/2:q/2
%         kk(:,2)=k2;
%         for k3=-q/2:q/2
%             kk(:,3)=k3;
%             tmp=-(sum((xx-kk).^2,2))/(4*b);
%             Q(:,k1+q/2+1,k2+q/2+1,k3+q/2+1)=exp(tmp)/c1;
%         end
%     end
% end

E=exp(b.*(2.*pi.*idx./(m*n)).^2);
E=E(:);

% The next block is faster
% E2=zeros(n,n,n);
% for k1=1:n    
%     for k2=1:n
%         for k3=1:n
%             E2(k1,k2,k3)=E(k1)*E(k2)*E(k3);
%         end
%     end
% end

EE=E*E.';
E3=zeros(n,n,n);
for k=1:n
    E3(:,:,k)=EE.*E(k);
end

u=beta.*E3;

low_idx_u=ceil((m*n-1)/2);
w=zeros(m*n,m*n,m*n);
w(low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1,...
  low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1,...
  low_idx_u+low_idx+1:low_idx_u+low_idx+1+n-1)=u;

W=fftshift(ifftn(ifftshift(w)));
W=W*numel(W);

% Slow version (was replaced by the subsequent block) 
% for k1=-q/2:q/2
%     for k2=-q/2:q/2
%         for k3=-q/2:q/2
%             idx=nu+repmat([k1 k2 k3],len_x,1);
%             idx=idx+repmat([low_idx_u low_idx_u low_idx_u],len_x,1);
%             idx=mod(idx,m*n);
%             for j=1:size(idx,1)
%                 g(j)=g(j)+Q(j,k1+q/2+1,k2+q/2+1,k3+q/2+1).*...
%                     W(idx(j,1)+1,idx(j,2)+1,idx(j,3)+1);
%             end
%         end
%     end
% end


c1=(2*sqrt(b*pi))^3;
xx=x.*m.*n/(2*pi)-nu;
offset1=zeros(len_x,3);
offset2=repmat([low_idx_u low_idx_u low_idx_u],len_x,1);
for k1=-q/2:q/2
    offset1(:,1)=k1;    
    for k2=-q/2:q/2
        offset1(:,2)=k2;
        for k3=-q/2:q/2
            offset1(:,3)=k3; %avoid allocating memory each time using repmat
            idx=nu+offset1;
            idx=idx+offset2;
            idx=mod(idx,m*n)+1;
            j=(idx(:,3)-1)*(n*m)^2 + (idx(:,2)-1)*n*m+idx(:,1); 
              % fast implementation of  sub2ind([m*n m*n m*n],idx(:,1),idx(:,2),idx(:,3));              
            tmp=-(sum((xx-offset1).^2,2))/(4*b);
            Q=exp(tmp);
            g=g+Q.*W(j);             
        end
    end
end
g=g./c1;