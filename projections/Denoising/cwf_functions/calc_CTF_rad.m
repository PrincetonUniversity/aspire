function [ctf_rad_all]=  calc_CTF_rad(use_CTF, n, index, ndef, def1,def2,B, lambda, r)

% Calculate CTF along a radial line (1D)
% resolution 2.82 same as add_CTF_env, this is arbitrary and should be
% changed to the CTF parameters used in add_CTF function
% Same alpha, other params
% Tejal Bhamre: updated Oct 2015

defocus=linspace(def1, def2, ndef);
ctf_rad_all=zeros(length(r),max(index));

if use_CTF==1
    for i=1:max(index)
        ctf_rad_all(:,i) = CTF_rad_cpyold(n, 2.82, lambda, defocus(:,i), 2, B, 0.07, r);
    end
    
elseif use_CTF==0
    ctf_rad_all=ones(length(r),max(index));
end
