function [pmin,fmin]=linmin(p,xi,f,data,tol,DEBUG)
%
% Given an n-dimensional point p and and n-dimensional direction xi, finds
% a point pmin where the function f takes a minimum along the direction xi
% from p. fmin is the value of f at pmin.
% data (optional) is additional info provided to the function f. Use [] to
% indicate no data.
%
% Yoel Shkolnisky, July 2008

if ~exist('DEBUG','var')
    DEBUG=0;
end

if ~exist('tol','var')
    tol=1.0e-8;
end

if ~exist('data','var')
    data=zeros(0);
end

[ax,bx,cx]=mnbrak(p,p+xi,f,data,0);
[pmin,fmin]=brentvec(ax,bx,cx,f,data,tol,DEBUG);

%% testing
%[xmin,fmin]=linmin(0,1,(@(x) x*(x-2)))
%[xmin,fmin]=linmin([0 0],[1 1],(@(x)cos(norm(x))),1.0e-8,1)
%[xmin,fmin]=linmin(-1,-0.9,(@(x) sin(10*x)),[],1.0e-8,1)
%[xmin,fmin]=linmin(0,1,(@(x) sin(0.5*x)+cos(x)),[],1.0e-8,1)
%[xmin,fmin]=linmin([1 -1],[0 1],@(x) 4.*(x(1))^2+(x(2))^2,[],1.e-8,1);


