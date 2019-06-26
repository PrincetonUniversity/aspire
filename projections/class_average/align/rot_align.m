function [corr, rot]=rot_align(m, coeff, list)
%   rotationally align pairs of images
%   Input:
%       m: frequencies computed from steerable basis
%       coeff: in cell structure, each cell contains coefficients with one
%       angular frequencies. 
%       list: of size Px2 where P is the number of pairs of images. 
%   Output:
%       rot: of size Px1 is the rotatinal alignment angle in degrees
%       ITER: the number of iterations
%   Zhizhen Zhao Feb 10 2012

N_theta=360;
P=size(list, 1);
max_m=max(m);
C=zeros(max_m, P);
m_list = 1:max_m;
m_list = m_list';
max_iter = 100;
precision = 1e-10;

for i=1:max_m+1
    C(i, :)=sum(conj(coeff{i}(:, list(1:P, 1))).*coeff{i}(:, list(1:P, 2)), 1);
end
C2=flipud(conj(C(2:end, :)));
B=real((2*max_m+1)*icfftd([C2; C], 1));
[~, rot]=max(B, [], 1);
rot=(rot-max_m-1)*N_theta/(2*max_m+1);

%%Newton-Raphson method for root finding
m_list_ang_1j = 1j * m_list*(pi/180);
constant_f_prime = repmat(m_list_ang_1j, 1, P).*C(2:end, :);
constant_f_prime2 = repmat(m_list_ang_1j .^ 2, 1, P).*C(2:end, :);
f_prime=@(x) sum(real(constant_f_prime.*exp(m_list_ang_1j * x)));
f_prime2=@(x) sum(real(constant_f_prime2.*exp(m_list_ang_1j * x)));

% Finding brackets, x1<x2 such that sign(f(x1)) != sign(f(x2)) and rot = (x1 + x2) / 2
step_size = 0.5;
x1 = rot;
x2 = rot;
bad_indices = ones(size(x1));
while any(bad_indices)
    x1(bad_indices) = x1(bad_indices) - step_size;
    x2(bad_indices) = x2(bad_indices) + step_size;
    f_x1 = f_prime(x1);
    f_x2 = f_prime(x2);
    bad_indices = f_x1 .* f_x2 > 0;
end

% Setting x1, x2 into x_low, x_high such that f(x_low)<f(x_high).
x_low = x1;
x_high = x2;
f_x_low = f_prime(x_low);
f_x_high = f_prime(x_high);
x_high_is_low = f_x_high < f_x_low;
tmp = x_low;
tmp(x_high_is_low) = x_high(x_high_is_low);
x_high(x_high_is_low) = x_low(x_high_is_low);
x_low = tmp;

% Handling f(x) = 0 case
f_x_low = f_prime(x_low);
f_x_low_0 = f_x_low == 0;
x_high(f_x_low_0) = x_low(f_x_low_0);
f_x_high = f_prime(x_high);
f_x_high_0 = f_x_high == 0;
x_low(f_x_high_0) = x_high(f_x_high_0);


rts = (x_low + x_high) / 2;
dx = abs(x_low - x_high);
dx_old = dx;
f = f_prime(rts);
df = f_prime2(rts);
for i=1:max_iter
    bisect_indices = ((rts - x_high) .* df - f) .* ((rts - x_low) .* df - f) > 0 | abs(2 * f) > abs(dx_old .* df);
    newton_indices = ~bisect_indices;
    dx_old = dx;

    % Handling out of range indices with Bisect step
    dx(bisect_indices) = (x_high(bisect_indices) - x_low(bisect_indices)) / 2;
    rts(bisect_indices) = x_low(bisect_indices) + dx(bisect_indices);

    % Handling the rest with newton step
    dx(newton_indices) = f(newton_indices) ./ df(newton_indices);
    rts(newton_indices) = rts(newton_indices) - dx(newton_indices);

    % Stop criteria
    if all(abs(dx) < precision)
        break
    end

    % Else update parameters
    f = f_prime(rts);
    df = f_prime2(rts);
    f_negative = f < 0;
    x_low(f_negative) = rts(f_negative);
    x_high(~f_negative) = rts(~f_negative);

    % Changing low and high of converged points
    converged = abs(dx) < precision;
    x_low(converged) = rts(converged);
    x_high(converged) = rts(converged);
end

rot=rts;
m_list = [0; m_list];
m_list_ang_1j = 1j * m_list*(pi/180);
C=C.*exp(m_list_ang_1j * rot);
corr=(real(C(1, :))+2*sum(real(C(2:end, :))))/2;

end
