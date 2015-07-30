function [xmin,fmin]=dbrent(ax,bx,cx,f,df,data,tol,DEBUG)
%
% Given a function f and its derivative df, and given a bracketing triplet
% ax, bx, cx (such that bx is between ax and cx, and f(bx) is less than
% bothf(ax) and f(cx)), the function isolates the minimum to a fractional
% precision tol using a modification of Brent's method that uses
% derivatives. The abscissa of the minimum is returned as xmin, and the
% value at the minumum is returned as fmin.
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

ITMAX=100;
ZEPS=1.0e-10;

if (numel(ax)~=1) || (numel(bx)~=1) || (numel(cx)~=1)
    error('ax, bx, and cx, muxt be scalars. Use dbrentvec instead');
end

a=min(ax,cx);
b=max(ax,cx);

a_sav=a;
b_sav=b;

e=0;
x=bx; w=bx; v=bx;
fx=myfeval(f,x,data); fv=fx; fw=fx;
dx=myfeval(df,x,data); dv=dx; dw=dx;

for iter=1:ITMAX
    xm=0.5*(a+b);
    tol1=tol*abs(x)+ZEPS;
    tol2=2.0*tol1;
    
    if (abs(x-xm) <= (tol2-0.5*(b-a))) 
        xmin=x;
        fmin=fx;
        
        if DEBUG
            displaymin(a_sav,b_sav,xmin,fmin,f,data);
        end
        
        return;
    end
    
    if (abs(e) > tol1)
        d1=2.0*(b-a);
        d2=d1;
        if (dw ~= dx) 
            d1=(w-x)*dx/(dx-dw); 
        end
        if (dv ~= dx) 
            d2=(v-x)*dx/(dx-dv); 
        end
        u1=x+d1;
        u2=x+d2;
        ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
        ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
        olde=e;
        e=d;
        if (ok1 || ok2)
            if (ok1 && ok2)
                if abs(d1) < abs(d2)
                    d=d1;
                else
                    d=d2;
                end                
            elseif (ok1)
                d=d1;
            else
                d=d2;
            end
            if (abs(d) <= abs(0.5*olde))
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=sign(tol1,xm-x);
                end
            else                
                if dx >= 0.0
                    e=a-x;
                else
                    e=b-x;
                end
                d=0.5*(e);
            end
        else            
            if dx >= 0.0
                e=a-x;
            else
                e=b-x;
            end
            d=0.5*(e);
        end
    else        
        if dx >= 0.0
            e=a-x;
        else
            e=b-x;
        end
        d=0.5*(e);
    end

    
    if (abs(d) >= tol1)
        u=x+d;
        fu=myfeval(f,u,data);
    else
        u=x+sign(tol1,d);
        fu=myfeval(f,u,data);
        if (fu > fx)
            xmin=x;
            fmin=fx;
            
            if DEBUG
                displaymin(a_sav,b_sav,xmin,fmin,f,data);
            end
            
            return;           
        end
    end
    
    du=myfeval(df,u,data);
    if (fu <= fx)
        if (u >= x) 
            a=x; 
        else 
            b=x; 
        end        
        v=w; fv=fw; dv=dw;    
        w=x;fw=fx; dw=dx;
        x=u; fx=fu; dx=du;
    else
        if (u < x) 
            a=u; 
        else 
            b=u; 
        end
        if (fu <= fw || w == x)
            v=w; fv=fw; dv=dw;
            w=u; fw=fu; dw=du;       
        elseif (fu < fv || v == x || v == w)
            v=u; fv=fu; dv=du;
        end
    end

end

warning('Too many iterations in dbrent');
fmin=0;
return;

%%
function  c=sign(a,b)
if b>=0
    c=abs(a);
else
    c=-abs(a);
end

%%
function displaymin(a,b,xmin,fmin,f,data)
Nt=200;
t=linspace(a,b,Nt);
vals=zeros(Nt,1);
for j=1:Nt
    vals(j)=myfeval(f,t(j),data);
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
%[xmin,fmin]=dbrent(ax,bx,cx,(@(x) x*(x-2)),(@(x) 2*x-2),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) cos(x),[],1);
%[xmin,fmin]=dbrent(ax,bx,cx,(@(x) cos(x)),(@(x) -sin(x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(10*x),[],1)
%[xmin,fmin]=dbrent(ax,bx,cx,(@(x) sin(10*x)),(@(x) 10*cos(10*x)),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(0.5*x)+cos(x),[],1)
%[xmin,fmin]=dbrent(ax,bx,cx,(@(x) sin(0.5*x)+cos(x)),(@(x) (0.5*cos(0.5*x)-sin(x))),[],1.0e-8, 1)

%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(x),[],1)
%[xmin,fmin]=dbrent(ax,bx,cx,(@(x) sin(x)), (@(x) (cos(x))),[],1.0e-8,1)

%[ax,bx,cx,fa,fb,fc,err]=mnbrak(0,1,@(x) besselj(0,x),[],1)
%[xmin,fmin]=dbrent(ax,bx,cx,(@(x) besselj(0,x)), (@(x) -besselj(1,x)),[],1.0e-8,1)
