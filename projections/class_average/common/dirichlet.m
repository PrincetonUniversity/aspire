% function y=dirichlet(t,m)
%
% compute the value of the dirichlet kernel of length m at point t
%
% Yoel Shkolnisky 22/10/01

function y=dirichlet(t,m)
y = zeros(size(t));
for k=1:prod(size(y))
    if (abs(t(k))<eps)
        y(k)=1;
    else
    y(k)= sin(pi*t(k))/(m*sin(pi*t(k)/m));
    end
end