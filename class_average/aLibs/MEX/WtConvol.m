function z=WtConvol(a,b,m)
% function z = WtConvol(a,b,m)
%
% A special weighted convolution function to accelerate
% the re-estimation step in HMM for motors.
% F. Sigworth 3 Apr 2007
% 
%  * Form the product c[w] = sum_u { a[u] * b[u+w] * m[u,u+w] }  
%  *  where the following sizes are assumed:
%  *        a is nu x 1,
%  *        b is nu x 1,
%  *        m is nu x nu.
% The convolution wraps around.
% 
% Typical use in the reestimation code is:
%     Xi=C.*WeightedConvol(alpha, beta, b);
% 
% Equivalent (slow) m-file code:
% 
nu=numel(a);
z=zeros(nu,1);
b=repmat(b,2,1); % duplicate to allow wrap-around
m=repmat(m,2,1);
for w=1:nu
    s=0;
    for u=1:nu
        s=s+a(u)*b(u+w-1)*m(u+w-1,u);
    end;
    z(w)=s;
end;
