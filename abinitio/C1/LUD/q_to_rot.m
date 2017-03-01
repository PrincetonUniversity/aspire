function rot = q_to_rot(a)

nrma = norm(a,'fro');
if abs(nrma - 1) > 1e-15;
    fprintf('norm(a)~=1');
    a = a./nrma;
end

rot = [a(1)^2+a(2)^2-a(3)^2-a(4)^2, 2*a(2)*a(3)-2*a(1)*a(4), 2*a(2)*a(4)+2*a(1)*a(3);
       2*a(2)*a(3)+2*a(1)*a(4), a(1)^2-a(2)^2+a(3)^2-a(4)^2, 2*a(3)*a(4)-2*a(1)*a(2);
       2*a(2)*a(4)-2*a(1)*a(3), 2*a(3)*a(4)+2*a(1)*a(2), a(1)^2-a(2)^2-a(3)^2+a(4)^2];
