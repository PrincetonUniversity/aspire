function [lambda,v,niter] = powerMethod(A)

[rows,cols] = size(A);
if rows ~= cols
    error('Input matrix must be square');
end

if A ~= transpose(A)
    error('Input matrix must be symmetric');
end

q      = randn(cols,1);
q_prev = zeros(cols,1);

niter = 0;
while norm(q-q_prev) > 1e-7    
    q_prev = q;
    q      = A*q_prev;
    q      = q/norm(q);
    niter = niter + 1;
end

lambda = q'*A*q;
v = q;