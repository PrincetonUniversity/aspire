function ciis_inds = compute_scls_inds(Ris_tilde,n_symm,n_theta)

nRisTilde = size(Ris_tilde,3);

% there are alltogether n_symm-1 self-cl in each image. If n_symm is even then the pairs of one of these
% lines (the one in the middle) is itself and is therefore useless (it is the one that lies on the global X_Y plane)
n_selfcl_pairs = floor((n_symm-1)/2); 

ciis = zeros(2*n_selfcl_pairs,nRisTilde,'uint16');
g = [cosd(360/n_symm) -sind(360/n_symm) 0; 
	 sind(360/n_symm)  cosd(360/n_symm) 0; 
	 0 				 0  1]; % a rotation of 360/n_symm degrees about the z-axis

for i=1:nRisTilde
        
    Ri_tilde = Ris_tilde(:,:,i);
    for s=1:n_selfcl_pairs
        Riigs = Ri_tilde.' * g^s * Ri_tilde;
        
        % extract a common-line index for each possible theta_ij
        c_1 = [-Riigs(8) ;  Riigs(7)]; %[-Riig(2,3)  Riig(1,3)];
        c_2 = [ Riigs(6) ; -Riigs(3)]; %[ Riig(3,2) -Riig(3,1)];
        
        cii1 = clAngles2Ind(c_1,n_theta);
        cii2 = clAngles2Ind(c_2,n_theta);
        
        if cii2 > n_theta/2
            cii2 = cii2 - n_theta/2;
            cii1 = cii1 + n_theta/2;
            cii1 = mod(cii1-1,n_theta)+1;
        end
        
        ciis(2*s-1,i)  = cii1;
        ciis(2*s,  i)  = cii2;
    end
end


ciis_inds = zeros(nRisTilde,n_selfcl_pairs);

for p = 1:n_selfcl_pairs
    ciis_inds(:,p) = sub2ind([n_theta,n_theta/2],ciis(2*p-1,:),ciis(2*p,:));
end

end