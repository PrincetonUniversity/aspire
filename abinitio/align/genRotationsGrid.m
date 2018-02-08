function [ rotations, angles ] = genRotationsGrid( resolution , JRJ, verbose)
%genRotationsGrid generate approximatly equally spaced rotations.
%   Input:
%       resolution - the number of samples per 2*pi.
%                    for example:
%                        resolution = 50  you get   4484 rotations
%                        resolution = 75  you get  15236 rotations
%                        resolution = 100 you get  39365 rotations
%                        resolution = 150 you get 129835 rotations
%       JRJ - (false by deafult) is wether to consider the 
%   Output:
%       rotations - 3X3Xnumber_of_rotations matrix. of all the rotations.
%       angles - 3Xnumber_of_rotations matrix. each column contains three
%                angles of the rotation in the following parametrization for
%                quaternions:
%                parametrization for SO3
%                   x = sin(tau)* sin(theta)* sin(phi);
%                   y = sin(tau)* sin(theta)* cos(phi);
%                   z = sin(tau)* cos(theta);
%                   w = cos(tau);
%
%   See also:
%            the function angs2q transfers angles into vector of quaternions
%
% Appriximated sampling error for given resolution:
%
% resolution | sampling error
%-----------------------------
%  10        | 1.495527e+00 
%  15        | 1.212038e+00 
%  20        | 7.507815e-01 
%  25        | 4.970845e-01 
%  30        | 5.207388e-01 
%  35        | 5.041540e-01 
%  40        | 3.398574e-01 
%  45        | 2.470507e-01 
%  50        | 3.010216e-01 
%  55        | 3.038057e-01 
%  60        | 1.934475e-01 
%  65        | 1.891904e-01 
%  70        | 2.041224e-01 
%  75        | 2.302468e-01 
%  80        | 1.624364e-01 
%  85        | 1.331552e-01 
%  90        | 1.567133e-01 
%  95        | 1.797956e-01 
% 100        | 1.085554e-01 


% for equaly spaced samples we need:
%    tau:   0 -> pi/2   (n/4 - points)
%    theta: 0 -> pi (n/2 * sin(tau) - points)
%    phi:   0 -> 2*pi (n*sin(tau)*sin(theta) - points)
%       
%
% in quaternions the JRJ duality is enterpretated as the the following
% equvalence relation:
% x     x
% y -> -y
% z -> -z
% w     w
%
% Thus to enumerate over all rotations up to JRJ one needs to enumerate
% only:
%    tau:   0 -> pi/2   (n/4 - points)
%    theta: 0 -> pi/2 (n/4 * sin(tau) - points)
%    phi:   0 -> 2*pi (n*sin(tau)*sin(theta) - points)
%

% Ver 1: Written in 2015.06.21 by Yariv Aizenbud
% Ver 2: added the angles variable and some documentation.
%        Written in 2015.06.22 by Yariv Aizenbud


if (nargin < 2)
    JRJ = false;
end

if (nargin < 3)
    verbose = false;
end
%tau_limit = pi/2;
%phi_limit = 2*pi;
% the code is duplicated. the only difference is on the for on theta
if JRJ == true   
    %theta_limit = pi/2;
    
    counter = 0;
    tau1_step = (pi/2)/(resolution/4);
    for tau1 = tau1_step/2:tau1_step:pi/2-tau1_step/2
        theta1_step = (pi/2)/(resolution/4*sin(tau1));
        for theta1 = theta1_step/2:theta1_step:pi/2-theta1_step/2
            phi1_step = (2*pi)/(resolution*sin(tau1)* sin(theta1));
            for phi1 = 0:phi1_step:2*pi - phi1_step
                counter = counter + 1;
            end
        end
    end
    if verbose
        log_message('There will be %d points', counter);
    end
    n_of_rotations = counter;


    % R will be all the possiable rotations as vectors of length 9

    if verbose
        log_message('Building rotations');
    end
    angles = zeros(3,n_of_rotations);
    rotations = zeros(3,3, n_of_rotations); 
    counter = 0;
    % possiable improvements:
    % 1. q_to_rot can be applied on a vector, maybe in GPU
    tau1_step = (pi/2)/(resolution/4);
    for tau1 = tau1_step/2:tau1_step:pi/2-tau1_step/2
        sintau1 = sin(tau1);
        costau1 = cos(tau1);
        theta1_step = pi/2/(resolution/4*sin(tau1));
        for theta1 = theta1_step/2:theta1_step:pi/2-theta1_step/2
            sintheta1 = sin(theta1);
            costheta1 = cos(theta1);
            phi1_step = (2*pi)/(resolution*sin(tau1)* sin(theta1));
            for phi1 = 0:phi1_step:2*pi - phi1_step
                counter = counter + 1;
                angles(:, counter) = [tau1, theta1, phi1];
                rotations(:,:,counter) = q_to_rot([sintau1* sintheta1* sin(phi1), sintau1* sintheta1* cos(phi1), sintau1* costheta1, costau1]);
            end
        end
    end
else % JRJ == flase
    %theta_limit = pi;

    counter = 0;
    tau1_step = (pi/2)/(resolution/4);
    for tau1 = tau1_step/2:tau1_step:pi/2-tau1_step/2
        theta1_step = (pi)/(resolution/2*sin(tau1));
        for theta1 = theta1_step/2:theta1_step:pi-theta1_step/2
            phi1_step = (2*pi)/(resolution*sin(tau1)* sin(theta1));
            for phi1 = 0:phi1_step:2*pi - phi1_step
                counter = counter + 1;
            end
        end
    end
    if verbose
        log_message('There will be %d points', counter);
    end

    n_of_rotations = counter;


    % R will be all the possiable rotations as vectors of length 9
    if verbose
        log_message('Building rotations');
    end
    angles = zeros(3,n_of_rotations);
    rotations = zeros(3,3, n_of_rotations); 
    counter = 0;
    % possiable improvements:
    % 1. q_to_rot can be applied on a vector, maybe in GPU
    tau1_step = (pi/2)/(resolution/4);
    for tau1 = tau1_step/2:tau1_step:pi/2-tau1_step/2
        sintau1 = sin(tau1);
        costau1 = cos(tau1);
        theta1_step = pi/(resolution/2*sin(tau1));
        for theta1 = theta1_step/2:theta1_step:pi-theta1_step/2
            sintheta1 = sin(theta1);
            costheta1 = cos(theta1);
            phi1_step = (2*pi)/(resolution*sin(tau1)* sin(theta1));
            for phi1 = 0:phi1_step:2*pi - phi1_step
                counter = counter + 1;
                angles(:, counter) = [tau1, theta1, phi1];
                rotations(:,:,counter) = q_to_rot([sintau1* sintheta1* sin(phi1), sintau1* sintheta1* cos(phi1), sintau1* costheta1, costau1].');
            end
        end
    end
end
