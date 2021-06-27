function err = evalO(R_true,R_est,G,X)

psi = X(1);
theta = X(2);
phi = X(3);
O = Rz(phi)*Ry(theta)*Rx(psi);

n = size(G,3);
dist = zeros(1,n);

for i = 1:n
    g = G(:,:,i);
    dist(1,i) = norm(R_true - O*g*O.'*R_est,'fro');
end
[err,~] = min(dist);

