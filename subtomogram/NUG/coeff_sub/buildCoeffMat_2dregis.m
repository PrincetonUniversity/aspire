%{
INPUT:
    data: data(:,:,i) is ith projection
    varNoise: variance of the noise
    bandlmt: bandlimit of the objective function
    t: truncate at (t+1)th degree representation
    numCompThreads: maximum number of CPU to be used
OUTPUT:
    coeffMat: Fourier coefficients of the objective function

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function coeffMat = buildCoeffMat_2dregis( sh , t )

const = 8*pi^2;
n = size(sh,2);
coeffMat = cell(t+1,1);

% when k = 0
Ck = zeros(n,n);
for i = 1:n
    for j = i+1:n
        Cij = sh(:,i)'*sh(:,i) + sh(:,j)'*sh(:,j) - 2*sh(1,i)*conj(sh(1,j))';
        Ck(i,j) = const*Cij;
    end
end
coeffMat{1} = Ck + ctranspose(Ck);

% other ks
for k = 1:t
    dk = 2*k + 1;
    Ck = zeros(n*dk,n*dk);
    idx = k^2; % number of sh coefficient in previous k
    for i = 1:n
        for j = i+1:n
            Cij = -2*sh(idx+1:idx+dk,i)*ctranspose(sh(idx+1:idx+dk,j));
            Ck(dk*(i-1)+1:dk*i,dk*(j-1)+1:dk*j) = const*Cij;
        end
    end
    coeffMat{k+1} = Ck + ctranspose(Ck);
end

fprintf('DONE!\n')

end









