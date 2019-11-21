function [ FBCoeff , B ] = extendFB( coeff_pos_k , B_pos )

ang_freq = B_pos(:,1);

FBCoeff = real(coeff_pos_k{1}); % FB coefficients
B = B_pos( ang_freq == 0 , : ); % Bessel stuff

for k = 1:max(ang_freq)
    
    % postive order
    FBCoeff = [ FBCoeff ; coeff_pos_k{k+1} ];
    
    B = [ B ; B_pos( ang_freq == k , : ) ];
    
    % negative order
    FBCoeff = [ FBCoeff ; conj(coeff_pos_k{k+1}) ];
    
    B_neg_k = B_pos( ang_freq == k , : );
    B_neg_k(:,1) = -B_neg_k(:,1);
    B = [ B ; B_neg_k ];
    
end

end









