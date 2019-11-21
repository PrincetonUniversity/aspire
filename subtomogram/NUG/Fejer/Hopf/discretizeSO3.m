function SO3 = discretizeSO3(t_S2)

load(['sphericalDesign/design',num2str(t_S2),'.mat']) % (t+1)^2 points
S2 = design;
S2_size = size(S2,1);

% discretize S^1
S1_size = round(sqrt(pi*S2_size));
alpha = 2*pi*linspace(0,1,S1_size); alpha = alpha(1:end);

% discretize S^2
[gamma,beta,~] = cart2sph(S2(:,1),S2(:,2),S2(:,3));
beta = pi/2 - beta; gamma = gamma + pi;

% SO(3) in Euler Z-Y-Z
SO3 = zeros( S2_size*S1_size , 3 );
count = 0;
for i = 1:S1_size
for j = 1:S2_size
    count = count + 1;
    SO3(count,:) = [ alpha(i) , beta(j) , gamma(j) ];
end
end

end









