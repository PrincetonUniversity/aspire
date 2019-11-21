function [ f ] = SquareCostF_subtomo(basis, sample_points, coeff_pos_k, ang_tilt, S, R12)
% compute the cost function for a pair of subtomogram slices with squared
% difference and Fourier Bessel representations.

% Input:
%   basis and sample_points: FB functions used and positon and weight in
%   [0, c].
%   coeff_pos_k: FB coefficients of two subtomograms, cell{i} stands for k=i-1
%   S: number of tilt series in one subtomogram
%   R12: the relative rotation between the i and j-th particle R2'*R1
% Output:
%   f: the scalar score between the two series of images
%
% Yuan Liu, 05/11/2016

% compute the relative Euler angles ZYZ between each pair of slices
% (a,b,g) = Rl'Rj'RiRk in a(k,l); 
a = zeros(S,S);
g = zeros(S,S);
R_tilt = zeros(3,3,S);
% compute from angles
for k = 1:S
    theta = (k - ceil(S/2))*ang_tilt;
    q = [cos(theta/2) sin(theta/2) 0 0]';
    R_tilt(:,:,k) = q_to_rot(q);
end
for k = 1:S
    for l = 1:S
        ang = eulang(R_tilt(:,:,l)'*R12*R_tilt(:,:,k)); % change/check the inverses?
        a(k,l) = ang(1);
        g(k,l) = ang(3);
    end
end

% compute cost coefficients b(k1,k2,q1,q2), where k starts from 0 and q
% starts from 1.
len = length(basis.ang_freqs);
f = 0;
for iter1 = 1:len
    k1 = basis.ang_freqs(iter1);
    q1 = basis.rad_freqs(iter1);
    for iter2 = 1:len  % can compute half of b since symmetric
        k2 = basis.ang_freqs(iter2);
        q2 = basis.rad_freqs(iter2);
        b = (sample_points.w .* basis.Phi_ns{k1+1}(:,q1))' * basis.Phi_ns{k2+1}(:,q2); 
        for pair = 1:S^2
            k = ceil(pair/S);
            l = mod(pair,S)+1;
            tik = -pi/2-g(k,l);
            tjl = pi/2+a(k,l);
            d1 = coeff_pos_k{k1+1}(q1,k)*exp(1i*k1*tik) ...
               - coeff_pos_k{k1+1}(q1,l+S)*exp(1i*k1*tjl);
            d2 = coeff_pos_k{k2+1}(q2,k)*exp(1i*k2*tik) ...
                - coeff_pos_k{k2+1}(q2,l+S)*exp(1i*k2*tjl);
            f = f+b*d1*conj(d2);
        end
    end
end


% radius of the disc containing the compact support of I
f = c*f;

end
