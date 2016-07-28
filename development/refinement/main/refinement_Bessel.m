function [corr, class, rot, shifts_b] = refinement_Bessel(data, ref, c, s, r_max )
%This is a function for refinement classification
%   Input: 
%       data: raw images
%       ref: references generated from the reconstructed model from last
%       iteration
%       c: CTF function
%       s: maximum shift range
%       r_max: the maximum radius of the particle
%   Output:
%       corr: normalized cross-correlation 
%       class: the class
% Zhizhen Zhao 2012

L=size(data, 1);
P_ref=size(ref, 3);
P=size(data, 3);

%add CTF to reference
for i=1:P_ref
    ref(:, :, i)=icfftn(c.*cfftn(ref(:, :, i)));
end;
ref=real(ref);
%%%%%%%%%%%

N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
ref=reshape(ref, L^2, P_ref);
ref=ref(r<=r_max, :);
[ b, R_freq ]=Bessel_ns_v5(r_max);
%% Including all shifted images
a=-s: 1: s; % checking shift for every other pixel;
num=length(a);
a1=repmat(a', num, 1);
a2=kron(a', ones(num, 1));
shifts=[a1, a2];
lshifts=size(shifts, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Change this part of the code to make it faster
coeff_ref_b = b\ref;
coeff_ref_b = coeff_ref_b(1:length(R_freq), :);
var_ref = var(coeff_ref_b, [], 2);
[~, id]=sort(var_ref, 'descend');
id=id(1:250);
coeff_ref_b = coeff_ref_b(id, :);
R_freq=R_freq(id);

b=b(:, id);
B=pinv(b);

coeff_b=zeros(size(b, 2), P, lshifts);

parfor k=1:lshifts
    data_shift= cryo_reshift_cart( data, -shifts(k, :) );
    data1=reshape(data_shift, L^2, P);
    data1=data1(r<=r_max, :);
    coeff = B*data1;
    coeff = coeff(1:length(R_freq), :);
    coeff_b(:, :, k)=coeff;
end;

clear data_shift data1

coeff_b=reshape(coeff_b, size(b, 2), P*lshifts);

coeff_b(R_freq==0, :)=coeff_b(R_freq==0, :)/sqrt(2);
coeff_ref_b(R_freq==0, :)=coeff_ref_b(R_freq==0, :)/sqrt(2);
for i=1:P*lshifts
    coeff_b(:, i)=coeff_b(:, i)/norm(coeff_b(:, i));
end;

for i=1:P_ref
    coeff_ref_b(:, i)=coeff_ref_b(:, i)/norm(coeff_ref_b(:, i));
end;
coeff_b(R_freq==0, :)=coeff_b(R_freq==0, :)*sqrt(2);
coeff_ref_b(R_freq==0, :)=coeff_ref_b(R_freq==0, :)*sqrt(2);

%Put Coeff in cell structure
k_max=max(R_freq);
Cell_coeff=cell(k_max+1, 1);
Cell_coeff_ref=cell(k_max+1, 1);
for i=1:k_max+1
    Cell_coeff{i}=coeff_b(R_freq==i-1, :);
end;

for i=1:k_max+1
   Cell_coeff_ref{i}=coeff_ref_b(R_freq==i-1, :);
end;

%Compute correlation, classification, alignment
[corr, rot, class, shifts_b ]=steerable_basis_corr_shift_v2( R_freq, Cell_coeff, Cell_coeff_ref, shifts );

end

