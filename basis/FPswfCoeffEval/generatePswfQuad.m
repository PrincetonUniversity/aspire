function [quadRulePtsX,quadRulePtsY,quadRuleWts,radialQuadPts,quadRuleRadialWts,numAngularPts] = generatePswfQuad(n,c,phi_approx_err,lambdaMax,epsilon,realFlag)
%% Generate radial quadrature nodes and weights
[radialQuadPts,radialQuadWts] = generatePswfRadialQuad(n,c,phi_approx_err,lambdaMax);

%% Find number of quadrature nodes per radius
% - Analytical solution
numAngularPts = ceil(exp(1)*radialQuadPts*c/2 + log(1/epsilon)) + 1;

% - Numerical solution, About 20% less points than the analytical solution above (for large c).
for i = 1:numel(radialQuadPts)
    angErrVec = abs(besselj(1:2*numAngularPts(i),c*radialQuadPts(i)));
    numAngularPts(i) = find( (sum(angErrVec)-cumsum(angErrVec)) < epsilon,1,'first');
    % - Make sure there is an even number of angular points for each radius
    % such that each node and its reflection are present, thus allowing for
    % x2 less computations for real-valued images.
    if (mod(numAngularPts(i),2)==1)
        numAngularPts(i) = numAngularPts(i) + 1;
    end
end

%% Generate Angular nodes and weights
quadRulePtsR = [];
quadRulePtsTheta = [];
quadRuleWts = [];
quadRuleRadialWts = zeros(numel(radialQuadPts),1);
for i=1:numel(radialQuadPts)        
    if (realFlag==1)
        quadRulePtsTheta = [quadRulePtsTheta 2*pi/numAngularPts(i)*(0:(numAngularPts(i)/2)-1)];
        quadRulePtsR = [quadRulePtsR repmat(radialQuadPts(i),1,numAngularPts(i)/2)];
    else
        quadRulePtsTheta = [quadRulePtsTheta 2*pi/numAngularPts(i)*(0:numAngularPts(i)-1)];
        quadRulePtsR = [quadRulePtsR repmat(radialQuadPts(i),1,numAngularPts(i))];
    end
    quadRuleRadialWts(i) = 2*pi/numAngularPts(i)*radialQuadPts(i)*radialQuadWts(i);
    quadRuleWts = [quadRuleWts repmat(quadRuleRadialWts(i),1,numAngularPts(i))];
end

quadRulePtsX = quadRulePtsR .* cos(quadRulePtsTheta);
quadRulePtsY = quadRulePtsR .* sin(quadRulePtsTheta);

end