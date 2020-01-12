
%Input: Grid is a 3xN array of points in the unit sphere
function [angular_dist]=getGridMeanDist(grid)

angDist=@(x,Y) (abs(acos(x*Y')))'*180/pi;
[~,D]=knnsearch(grid',grid','Distance',angDist,'k',4);
D=D(:,2:end);
angular_dist=mean(D(:));

end