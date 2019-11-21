function [ scr ] = cl_to_score( ind, corrstack, n_theta )
% find the corresponding cross-correlation of a common line configuration
%
% ind is a 2ns*2ns matrix representing the common line numbers of the
% column number slice on the row number slice, where the numbers have an
% orientation looking from the -x axis.
%
% corrstack is the output of commonlines_gaussian_vsub.m where the first
% dim is the cross-correlation table with linear subscription of length
% n_theta * (n_theta/2).
%
% n_theta is the angular resolution.

ns = length(ind)/2;
scr = 0;
for i = 1:ns % index of the first subtomogram
    for j = 1:ns
        I = ind(i,j+ns);
        J = ind(j+ns,i);
        if J > n_theta/2
            J = J - n_theta/2;
            I = I + n_theta/2*sign(n_theta/2-I);
        end
        line = sub2ind([n_theta,n_theta/2], I, J);
        s = corrstack(line,i,j+ns);
        scr = scr + s;
    end
end
end

