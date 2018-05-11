function [lower_index,upper_index] = bsearch(x,LowerBound,UpperBound)
% BSEARCH Binary search in a sorted vector.
%
% Binary O(log2(N)) search of the range of indices of all elements of x
% between LowerBound and UpperBound. If no elements between LowerBound and
% Upperbound are found, the returned lower_index and upper_index are empty.
% The array x is assumed to be sorted from low to high, and is NOT verified
% for such sorting.
%
% Input parameters:
% x             A vector of sorted values from low to high.       
% LowerBound    Lower boundary on the values of x in the search.
% UpperBound    Upper boundary on the values of x in the search.
%
% Output parameters:
% lower_index   The smallest index such that
%               LowerBound<=x(index)<=UpperBound 
% upper_index   The largest index such that
%               LowerBound<=x(index)<=UpperBound
%
% Examples:
%   [startidx,endidx]=bsearch([1 1.5 2 2 2 3],2,2)
%   [startidx,endidx]=bsearch([1 1.5 2 2 2 3],1.5,2)
%
% Based on code from 
% http://stackoverflow.com/questions/20166847/faster-version-of-find-for-sorted-vectors-matlab
%
% Yoel Shkolnisy, October 2014.

if LowerBound>x(end) || UpperBound<x(1) || UpperBound<LowerBound
    % no indices satify bounding conditions
    lower_index = [];
    upper_index = [];
    return;
end

lower_index_a=1;
lower_index_b=length(x); % x(lower_index_b) will always satisfy lowerbound
upper_index_a=1;         % x(upper_index_a) will always satisfy upperbound
upper_index_b=length(x);

%
% The following loop increases _a and decreases _b until they differ 
% by at most 1. Because one of these index variables always satisfies the 
% appropriate bound, this means the loop will terminate with either 
% lower_index_a or lower_index_b having the minimum possible index that 
% satifies the lower bound, and either upper_index_a or upper_index_b 
% having the largest possible index that satisfies the upper bound. 
%
while (lower_index_a+1<lower_index_b) || (upper_index_a+1<upper_index_b)

    lw=floor((lower_index_a+lower_index_b)/2); % split the upper index

    if x(lw) >= LowerBound
        lower_index_b=lw; % decrease lower_index_b (whose x value remains \geq to lower bound)   
    else
        lower_index_a=lw; % increase lower_index_a (whose x value remains less than lower bound)
        if (lw>upper_index_a) && (lw<upper_index_b)
            upper_index_a=lw;% increase upper_index_a (whose x value remains less than lower bound and thus upper bound)
        end
    end

    up=ceil((upper_index_a+upper_index_b)/2);% split the lower index
    if x(up) <= UpperBound
        upper_index_a=up; % increase upper_index_a (whose x value remains \leq to upper bound) 
    else
        upper_index_b=up; % decrease upper_index_b
        if (up<lower_index_b) && (up>lower_index_a)
            lower_index_b=up;%decrease lower_index_b (whose x value remains greater than upper bound and thus lower bound)
        end
    end
end

if x(lower_index_a)>=LowerBound
    lower_index = lower_index_a;
else
    lower_index = lower_index_b;
end
if x(upper_index_b)<=UpperBound
    upper_index = upper_index_b;
else
    upper_index = upper_index_a;
end

if upper_index<lower_index
    % The requested range of keys not found in the array.
    % Return empty lower and uppoer bounds.
    lower_index = [];
    upper_index = [];
end
    