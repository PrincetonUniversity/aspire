function q=qrand(K)
%
% q=qrand(K)
%
% Generate K random uniformly distributed quaternions.
% Each quaternions is a four-elements column vector. Returns a matrix of
% size 4xK.
%
% The 3-sphere S^3 in R^4 is a double cover of the rotation group SO(3),
% SO(3) = RP^3. 
% We identify unit norm quaternions a^2+b^2+c^2+d^2=1 with group elements. 
% The antipodal points (-a,-b,-c,-d) and (a,b,c,d) are identified as the
% same group elements, so we take a>=0.
%

q = randn(4,K);
%q(:,1) = [0;1;0;1];
l2_norm = sqrt(q(1,:).^2 + q(2,:).^2 + q(3,:).^2 + q(4,:).^2);
for i=1:4
    q(i,:) = q(i,:) ./ l2_norm;
end;
for k=1:K
    if (q(1,k) < 0)
        q(:,k) = -q(:,k);
    end;
end;
