function [corr, rot, class, shift, ITER]=steerable_basis_corr_shift_v2(m, coeff, coeff_ref, shifts)

N_theta=360;
lshifts=length(shifts);
P=size(coeff{1}, 2);
P=P/lshifts;
P_ref=size(coeff_ref{1}, 2);
max_m=max(m);
rot=zeros(P, 1);
class=zeros(P, 1);
corr=zeros(P, 1);
m_list=[1:max_m]';
C=zeros(max_m+1, P_ref*lshifts);  
a=sqrt(lshifts);
shift=zeros(P, 2);
ITER=zeros(P, 1);

ps=matlabpool('size');
if ps==0
    matlabpool open
end
parfor k=1:P  
    [ C ] = make_C(coeff, coeff_ref, k, max_m, lshifts);
    C2=flipud(conj(C(2:end, :)));
    B=real((2*max_m+1)*icfft_1dstack([C2; C]));
    [corr(k), tmp]=max(abs(B(:)));
    [rot_tmp, id]=ind2sub([2*max_m+1, P_ref*lshifts], tmp);
    [id2, class(k)]=ind2sub([lshifts, P_ref], id);
    
    shift(k, :)=shifts(id2, :); 
    rot_tmp=(rot_tmp-max_m-1)*N_theta/(2*max_m+1);
    
    %%%newton method for root finding 
    x_old=-1;
    x_new=rot_tmp;
    precision=0.001;
    f_prime=@(x) sum(real(sqrt(-1)*m_list*(pi/180).*C(2:end, id).*exp(sqrt(-1)*m_list*x*pi/180)));
    f_prime2=@(x) sum(real((-1)*m_list.^2*(pi/180)^2.*C(2:end, id).*exp(sqrt(-1)*m_list*x*pi/180)));
    iter=0;
    while abs(x_new - x_old) > precision
       x_old = x_new;
       delta=f_prime(x_old)/f_prime2(x_old);
       if abs(delta)>10
           delta=sign(delta)*10*rand(1);
       end;
       x_new = x_old - delta;
       iter=iter+1;
       if iter>10, break, end
    end;
    ITER(k)=iter;
    rot(k)=x_new;
end

matlabpool close;

end

function [ C ] = make_C(coeff, coeff_ref, k, max_m, lshifts)
P_ref=size(coeff_ref{1}, 2);
P=size(coeff{1}, 2)/lshifts;
C=zeros(max_m+1, P_ref*lshifts);  
for i=1:max_m+1
    tmp = coeff{i}(:, [k:P: P*(lshifts-1)+k])'*coeff_ref{i};
    C(i, :)=reshape(tmp, 1, P_ref*lshifts);
end;

end
