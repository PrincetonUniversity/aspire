function [pmin,fmin,iter]=frprmn(p,f,df,data,ftol,DEBUG)
%
% Given a starting point p, run Fletcher-Reeves-Polak-Ribirer minimization
% on the function f using its gradient as calculated by df. The convergence
% tolerance in ftol. Returns the location of the minimum pmin, the value at
% the minimum fmin, and the number of iterations performed iter.
%
% Yoel Shkolnisky, July 2008.

if ~exist('DEBUG','var')
    DEBUG=0;
end

if ~exist('ftol','var')
    ftol=1.0e-8;
end

if ~exist('data','var')
    data=zeros(0);
end


ITMAX=200;
EPS=1.0e-10;
TOL=2.0e-8;  % Square root of machine precision

fp=myfeval(f,p,data);
xi=myfeval(df,p,data);
g=-xi;
h=g; 
xi=h;

PLOT_PROGRESS=0;
if DEBUG
    if prod(size(p))==2
        PLOT_PROGRESS=1;
        fig_h=gcf;
        clf;
        subplot(1,2,1);
        plotmap(p,[0 0],f,data,0);
        hold on;
        scatter(p(1),p(2),30,'r','filled');
        title('local map');
        hold off;
        
        subplot(1,2,2);
        plotmap(p,[0 0],f,data,1);
        hold on;
        scatter(p(1),p(2),30,'m','filled');
        title('prgress');
        hold off;
        
        p_old=p;
    end
end

for its=1:ITMAX
    iter=its;   
    [p,fret]=dlinmin(p,xi,f,df,data,TOL,0);
    %[p,fret]=linmin(p,xi,f,data,TOL,0);
    
    if PLOT_PROGRESS
        figure(fig_h);
        subplot(1,2,1);
        plotmap(p,p_old,f,data,0);
        hold on;
        scatter(p(1),p(2),30,'r','filled');
        title('local map');
        hold off;
        
        subplot(1,2,2);
        hold on;
        line([p_old(1) p(1)],[p_old(2) p(2)],'LineStyle',':');
        p_old=p;
        scatter(p(1),p(2),15,'r','filled');
        hold off;

    end
    
    if (2*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS))
        pmin=p;
        fmin=fret;
        return
    end
    
    fp=fret;
    xi=myfeval(df,p,data);
    gg=sum(g.^2);
    dgg=sum((xi+g).*xi);
    
    if (gg==0)
        pmin=p;
        fmin=fret;
        return
    end
    
    gam=dgg/gg;    
    g=-xi;
    h=g+gam*h;
    xi=h;    
end

warning('Too many iterations in frprmn');
pmin=p;
fmin=fret;

function plotmap(p,p_old,f,data,fixedview)
if fixedview    
    [X,Y]=meshgrid(linspace(-10,10,100),linspace(-5,5,100));
else
    dx=5*abs(p(1)-p_old(1)); dy=2*abs(p(2)-p_old(2));
    d=max(dx,dy);
    [X,Y]=meshgrid(p(1)+linspace(-d,d,100),p(2)+linspace(-d,d,100));    
end
vals=zeros(size(X));
for j=1:length(X(:))
    vals(j)=myfeval(f,[X(j) ; Y(j)],data);
end
contourf(X,Y,vals,50);
colormap(cool);


%% testing
%[xmin,fmin]=frprmn([1 2],(@(x) x(1).^2+x(2).^2), (@(x) [2*x(1),2*x(2)] ),[],1.0e-8, 1)
%[xmin,fmin]=frprmn([0 10],(@(x) 100*x(1).^2+x(2).^2), (@(x) [200*x(1),2*x(2)] ),[],1.0e-8,0)
%[xmin,fmin]=frprmn([1 1.5],(@(x) 3*x(1).^2+x(2).^2), (@(x) [6*x(1),2*x(2)] ),[],1.0e-8,1)
%[xmin,fmin]=frprmn([1 1],(@(x) x(1)^2+(x(1)^2+x(2))^2), (@(x) [2*x(1)+4*x(1)*(x(1)^2+x(2)), 2*(x(1)^2+x(2))] ),[],1.0e-8,1)
%[xmin,fmin]=frprmn([4 3],(@(x) (x(2)-sin(10*x(1)))^2 + x(1)^2+x(2)^2), (@(x) [-20*cos(10*x(1))*(x(2)-sin(10*x(1)))+2*x(1), 2*(x(2)-sin(10*x(1)))+2*x(2)]),[],1.0e-8,1)
