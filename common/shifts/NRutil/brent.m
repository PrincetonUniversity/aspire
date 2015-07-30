function [xmin,fmin]=brent(ax,bx,cx,f,data,tol,DEBUG)
%
% Given a function handle f, and given a bracketing triplet of abscissas
% ax, bx, cx (sych that bx is between ax and cx), and f(bx) is less than
% both f(ax) and f(cx), this function isolates the minimum to a fractional
% precision of about tol using Brent's method. The abscissa of the minimum
% is returned as xmin, and the minimum function value is returned as fmin.
% data (optional) is additional info provided to the function f. Use [] to
% indicate no data.
%
% Yoel Shkolnisky, July 2008.

if ~exist('DEBUG','var')
    DEBUG=0;
end

if ~exist('tol','var')
    tol=1.0e-8;
end

if ~exist('data','var')
    data=zeros(0);
end

ITMAX=100;            % Maximum allowed number of iterations.
CGOLD=(3-sqrt(5))/2;  % Golden ratio.
ZEPS=1.0e-10;         % Protects against trying to achieve fractional 
                      % accuracy for a minimum that is exactly zero.

                      
if (prod(size(ax))~=1) || (prod(size(bx))~=1) || (prod(size(cx))~=1)
    error('ax, bx, and cx, muxt be scalars. Use brentvec instead');
end
                      
e=0.0;

a=min(ax,cx);
b=max(ax,cx);

a_sav=a;
b_sav=b;

x=bx; w=bx; v=bx;
fx=myfeval(f,x,data); fv=fx; fw=fx;

for iter=1:ITMAX
    xm=0.5*(a+b);
    tol1=tol*abs(x)+ZEPS;
    tol2=2.0*(tol1);
    if (abs(x-xm) <= (tol2-0.5*(b-a)))
        xmin=x;
        fmin=fx;                
        
        if DEBUG
            displaymin(a_sav,b_sav,xmin,fmin,f,data);
        end

        return;
    end
    
    if (abs(e) > tol1)
        r=(x-w)*(fx-fv);
        q=(x-v)*(fx-fw);
        p=(x-v)*q-(x-w)*r;
        q=2.0*(q-r);
        if (q > 0.0)
            p = -p;
        end
        q=abs(q);
        etemp=e;
        e=d;
        if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            if x>=xm
                e=a-x;
            else
                e=b-x;
            end
            d=CGOLD*(e);
        else
            d=p/q;
            u=x+d;
            if (u-a < tol2 || b-u < tol2)
                d=sign(tol1,xm-x);
            end
        end
    else
        if x >= xm
            e=a-x;
        else
            e=b-x;
        end
        d=CGOLD*(e);
    end

    if abs(d) >= tol1
        u=x+d;
    else
        u=x+sign(tol1,d);
    end

    fu=myfeval(f,u,data);
    if (fu <= fx)
        if (u >= x) 
            a=x; 
        else
            b=x; 
        end
        v=w; w=x; x=u;
        fv=fw; fw=fx; fx=fu;
    else
        if (u < x) 
            a=u; 
        else
            b=u; 
        end;
        if (fu <= fw || w == x)
            v=w;
            w=u;
            fv=fw;
            fw=fu;
        elseif (fu <= fv || v == x || v == w)
            v=u;
            fv=fu;
        end
    end  
end

warning('Too many iterations in brent');
xmin=x;
fmin=fx;

if DEBUG
    displaymin(a_sav,b_sav,xmin,fmin,f,data);
end

return;


function  c=sign(a,b)
if b>=0
    c=abs(a);
else
    c=-abs(a);
end


function displaymin(a,b,xmin,fmin,f,data)
Nt=200;
t=linspace(a,b,Nt);
vals=zeros(Nt,1);
for j=1:Nt
    vals(j)=feval(f,t(j));
end

fa=myfeval(f,a,data);
fb=myfeval(f,b,data);

plot(t,vals);
hold on;
scatter(a,fa,15,'r','filled');
scatter(b,fb,15,'r','filled');
scatter(xmin,fmin,15,'r','filled');
text(a,fa+0.1,'a');
text(b,fb+0.1,'b');
text(xmin,fmin+0.1,'x');
hold off;


%% testing
%[ax,bx,cx,fa,fb,fc,err]=mnbrak(0.1,3,@(x) x*(x-2),[],1);
%[xmin,fmin]=brent(ax,bx,cx,(@(x) x*(x-2)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) cos(x),[],1);
%[xmin,fmin]=brent(ax,bx,cx,(@(x)cos(x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(10*x),[],1)
%[xmin,fmin]=brent(ax,bx,cx,(@(x) sin(10*x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(0.5*x)+cos(x),[],1)
%[xmin,fmin]=brent(ax,bx,cx,(@(x) sin(0.5*x)+cos(x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(x),[],1)
%[xmin,fmin]=brent(ax,bx,cx,(@(x) sin(x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc,err]=mnbrak(0,1,@(x) besselj(0,x),[],1)
%[xmin,fmin]=brent(ax,bx,cx,(@(x) besselj(0,x)),[],1.0e-8,1)
