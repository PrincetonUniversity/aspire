function [c,ceq] = Bfitconstraints(params,xData,yData)
c = yData - params(1).*exp(-params(2).*xData.^2);
ceq = [];
