function s_Wishart = KN_s_Wishart(alpha,beta)
% function s_Wishart = KN_s_Wishart(alpha,beta)
%
% Code by Shira Kritchman and Boaz Nadler
% 2008, Weizmann Institute of Science
% --------------------------------------------
% DESCRIPTION:
%   This function computes an approximate inverse of the 
%   TW (Tracy-Widom) distribution F_beta: 
%   if X ~ F_beta, then the function returns a value s_Wishart such that 
%   Pr{X > s_Wishart} ~ alpha/100.
%   This value is used in the algorithm for rank estimation, KN_rankEst.
%
% INPUT:
%   alpha   -  confidence level (given in percentage)
%   beta    -  indicator for real (1) or complex (2) TW distribution
%
% OUTPUT:
%   s_Wishart - threshold value for TW distribution with confidence level alpha
% --------------------------------------------
% FOR MORE DETAILS SEE:
%   S. Kritchman and B. Nadler, Determining the number of components in a factor model
%   from limited noisy data, 2008
%   ------------------------------------------
% Here we use the asymptotics of the TW distribution for large values of x
% These can be found, for example, in 
% J. Baik, R. Buckingham and J. DiFranco,
% Asymptotics of the Tracy-Widom distributions and the total integral of a Painleve II function,
% Comm. Math. Phys., vol. 280, no. 2, pp. 463--497, 2008.
% --------------------------------------------

if beta == 2
    s_Wishart = (-3/4 * log(16 * pi * alpha/100) )^(2/3); 
%    fprintf('s_Wishart %f %f\n',s_Wishart,alpha);  pause; 
    return; 
end

if beta == 1
    s_Wishart = (-3/2 * log(4*sqrt(pi) * alpha/100 ))^(2/3); 
    return; 
end

