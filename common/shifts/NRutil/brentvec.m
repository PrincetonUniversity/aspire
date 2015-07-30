function [xmin,fmin]=brentvec(ax,bx,cx,f,data,tol,DEBUG)
%
% Like brent but ax, bx, and cx are colinear vectors, and f is a scalar
% funtion on R^{n}.
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
EPS=1.0e-14;                      

if isa(ax,'single') || isa(bx,'single') || isa(cx,'single')   
    ZEPS=1.0e-3;
    EPS=1.0e-6;
end

x0=ax;
dvec=(bx-ax);
if norm(dvec)<EPS % no direction is given. The minimum is the current point
    xmin=x0;
    fmin=myfeval(f,xmin,data);
    return;
end
   
dvec=dvec/norm(dvec);
if (abs(dot(cx-ax,dvec))/(norm(dvec)*norm(cx-ax))-1) > EPS
    error('ax, bx, cx must be colinear');
end

% convert ax, bx, and cx to time coordinates on the 1D coordinate system
% with origin at x0 and positive direction dvec.
ax=0;
bx=dot(bx-x0,dvec)/norm(dvec).^2;
cx=dot(cx-x0,dvec)/norm(dvec).^2;

a=min(ax,cx);
b=max(ax,cx);

a_sav=a;
b_sav=b;

e=0.0;
x=bx; w=bx; v=bx;
fx=myfeval(f,point(x0,dvec,x),data); fv=fx; fw=fx;

for iter=1:ITMAX
    xm=0.5*(a+b);
    tol1=tol*abs(x)+ZEPS;
    tol2=2.0*(tol1);
    if (abs(x-xm) <= (tol2-0.5*(b-a)))
        xmin=point(x0,dvec,x);
        fmin=fx;
        
        if DEBUG
            displaymin(x0,dvec,a_sav,b_sav,xmin,fmin,f,data);
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

    fu=myfeval(f,point(x0,dvec,u),data);
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
xmin=point(x0,dvec,x);
fmin=fx;

if DEBUG
    displaymin(x0,dvec,a_sav,b_sav,xmin,fmin,f,data);
end

return;


function  c=sign(a,b)
if b>=0
    c=abs(a);
else
    c=-abs(a);
end

function x=point(x0,dvec,t)
x=x0+t.*dvec;

function displaymin(x0,dvec,a,b,xmin,fmin,f,data)
Nt=200;
t=linspace(a,b,Nt);
vals=zeros(Nt,1);
for j=1:Nt
    vals(j)=myfeval(f,point(x0,dvec,t(j)),data);
end

fa=myfeval(f,point(x0,dvec,a),data);
fb=myfeval(f,point(x0,dvec,b),data);

xmin=dot(xmin-x0,dvec)./norm(dvec).^2;

plot(t,vals);
hold on;
scatter(a,fa,15,'r','filled');
scatter(b,fb,15,'r','filled');
scatter(xmin,fmin,15,'r','filled');
text(double(a),double(fa)+0.1,'a'); % Text has bug and crashes if the input is single.
text(double(b),double(fb)+0.1,'b');
text(double(xmin),double(fmin)+0.1,'x');
hold off;


%% testing
%[ax,bx,cx,fa,fb,fc,err]=mnbrak(0.1,3,@(x) x*(x-2));
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x) x*(x-2)),[],1.0e-8)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) cos(x));
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x)cos(x)))

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(10*x))
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x) sin(10*x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(0.5*x)+cos(x))
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x) sin(0.5*x)+cos(x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(x))
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x) sin(x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc,err]=mnbrak(0,1,@(x) besselj(0,x))
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x) besselj(0,x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc,err]=mnbrak([1 0 0],[1 1 1],@(x) norm(x)*norm(x-2))
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x) norm(x)*(norm(x-2))),[],1.0e-8,1)


%[ax,bx,cx,fa,fb,fc,err]=mnbrak([0 0 0],[1 1 0],@(x) besselj(0,norm(x)))
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x) besselj(0,norm(x))),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak([0 0],[1 2],@(x) cos(norm(x)));
%[xmin,fmin]=brentvec(ax,bx,cx,(@(x)cos(norm(x))),[],1.0e-8,1)
