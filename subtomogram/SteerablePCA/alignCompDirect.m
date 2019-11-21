function [eul1 ] = alignCompDirect( vecSt, vecS, siz_n, bandl )
% Align two eigenvectors representing the principal component in expansion
% coefficient.

% Input
% vec1: coefficient to align
% vec2: coefficient for reference volume

% Output
% rotated_vec1: rotated coefficient of vec1 according to vec2
N = 2*bandl-1;
alpha = pi*(0:N)/bandl;
beta = pi*((0:N)*2+1)/(4*bandl);
[alpha,beta,gamma] = ndgrid(alpha,beta,alpha);
angles = [alpha(:) beta(:) gamma(:)];
A = vec2cell(vecSt,siz_n);
AN = cell(size(siz_n,1), 1);
ANvec = zeros(size(vecSt,1),(2*bandl)^3);
maxL = size(siz_n,1)-1;
for i = 1:(2*bandl)^3
    angle = angles(i,:)';
    for ll = 0:maxL
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        %[T,Tinv] = realY_to_complexY(ll);
        %Wl = real(Tinv*Wl*T);
        al = A{ll+1,1}*conj(Wl); % give complex numbers, use real Y?
        AN{ll+1,1} = al;
    end
    ANvec(:,i) = cell2vec(AN);
end
diff = zeros((2*bandl)^3,1);
for i = 1:(2*bandl)^3
    diff(i) = norm(ANvec(:,i)-vecS);
end
[~,ind] = min(diff);
eul1 = angles(ind,:);

end

