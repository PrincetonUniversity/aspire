%{
INPUT:
    w: weights for numerical integration
    P: projections in Fourier space
    theta: angle discretization
    alpha: 1st Euler angle
    gamma: 3rd Euler angle
OUTPUT:
    objVal: value of the objective function at (alpha,gamma)X(i,j)

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function objVal = getObjVal_inBatches( w , P , idx_I , idx_J , objVal_subI , objVal_subJ )

% L1 norm
objVal = abs( P( : , idx_I , objVal_subI ) - P( : , idx_J , objVal_subJ ) );
objVal = sum( objVal .* repmat( w , 1 , size(objVal,2) , size(objVal,3) ) , 1 );
objVal = squeeze(objVal);

end









