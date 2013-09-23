function g=nufft_t_3d_execute_2(x,prepdata)
%
% See nufft_t_3d_prepare_2.
%
% Yoel Shkolnisky, February 2010.

g=nufftt3dexecutemx_2(x,prepdata);

%
% This is a slower implementation of the above MEX file:
%
% if length(size(x))~=2
%     error('x must be a mx3 array')
% end
% 
% if size(x,2)~=3
%     error('x must be a mx3 array');
% end
% 
% 
% W=prepdata.W;
% b=prepdata.b;
% m=prepdata.m;
% q=prepdata.q;
% n=prepdata.n;
% 
% len_x=size(x,1);
% g=zeros(len_x,1);
% nu=round(x*m*n/(2*pi));
% low_idx_u=ceil((m*n-1)/2);
% 
% c1=(2*sqrt(b*pi))^3;
% xx=x.*m.*n/(2*pi)-nu;
% offset1=zeros(len_x,3);
% offset2=repmat([low_idx_u low_idx_u low_idx_u],len_x,1);
% for k1=-q/2:q/2
%     offset1(:,1)=k1;    
%     for k2=-q/2:q/2
%         offset1(:,2)=k2;
%         for k3=-q/2:q/2
%             offset1(:,3)=k3; %avoid allocating memory each time using repmat
%             idx=nu+offset1;
%             idx=idx+offset2;
%             idx=mod(idx,m*n)+1;
%             j=(idx(:,3)-1)*(n*m)^2 + (idx(:,2)-1)*n*m+idx(:,1); 
%               % fast implementation of  sub2ind([m*n m*n m*n],idx(:,1),idx(:,2),idx(:,3));
%               
%              tmp=-(sum((xx-offset1).^2,2))/(4*b);
%              Q=exp(tmp);
%              g=g+Q.*W(j);             
%         end
%     end
% end
% g=g./c1;
