function [data]=refinement_align(rot, shifts_b, data)

% L=size(data, 1);
% step_size=1;
P=size(data, 3);
% H=sparse(1:P, 1:P, exp(sqrt(-1)*rot*pi/180), P, P);
% s_tmp=step_size*(shifts_b(:, 1)+sqrt(-1)*shifts_b(:, 2));
% s_2=H*s_tmp;
% 
% shifts_b(:, 1)=real(s_2);
% shifts_b(:, 2)=imag(s_2);

data=shift_images(data, shifts_b);
  
parfor k=1:P
    data(:, :, k)=fastrotate(data(:, :, k), -rot(k));
end;
