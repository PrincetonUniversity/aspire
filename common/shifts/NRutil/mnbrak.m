function [ax,bx,cx,fa,fb,fc,err]=mnbrak(ax,bx,f,data,DEBUG)
%
% Bracket a minimum of the function f using the initial abscissas ax and
% bx.
% The function searches in the downhill direction and returns three points
% that bracket a minimum of the function. The function also returns the
% values of f at these points. f is a handle to a scalar function that take
% an n-dimenional point.
%
% ax and bx are the initial search points, each in dimension n. The
% function performs a one dimensional search on the line between ax and bx.
%
% err=0 if everything went OK. err=1 if the bracketing points become to far
% apart. This probably means that the function is deacresing.
%
% ax, bx, cx are the bracketing points. fa, fb, fc are the corresponding
% function values.
%
% data (optional) is additional info provided to the function f. Use [] to
% indicate no data.
%
% Yoel Shkolnisky, July 2008.

if ~exist('DEBUG','var')
    DEBUG=0;
end

if ~exist('data','var')
    data=zeros(0);
end

GOLD=2/(3-sqrt(5))-1;
GLIMIT=100;
TINY=1.0e-20;
TLIMIT=1e6;
n=length(ax);

if length(bx)~=n
    error('ax and bx must have same dimensions');
end
err=0;

ta=0;     % Minimize on the line between ax and bx. ax is considered 0 (v0),
tb=norm(bx-ax);     % bx is considered 1, and the vector bx-ax is one step.

if tb==0 % ax and bx are the same. No direction is given.
    cx=bx;
    return;
end

v0=ax;    % Once one-dimensional bracketing is determined, the corresponding 
vec=bx-ax;  % points in n-dimensional space are computed.
vec=vec/norm(vec);

fa=myfeval(f,ax,data);
fb=myfeval(f,bx,data);

if (fb>fa) % Switch ta and tb so that we can go downhill in the direction from ta to tb.
    tmp=ta; ta=tb; tb=tmp;
    tmp=fa; fa=fb; fb=tmp;
end

tc=tb+GOLD*(tb-ta); % First guess for tc.
fc=myfeval(f,v0+tc*vec,data);

done=(fc >= fb);
while (~done)    %Keep returning until we bracket. 
    
    if abs(tc)>TLIMIT  % don't allow the farthest point to go too far
        err=1;
        done=1;
        continue;
    end
    
    r=(tb-ta).*(fb-fc);   % Compute u by parabolic extrapolation from 
    q=(tb-tc).*(fb-fa);   % ta, tb, and tc. Tiny is used to prevent division by zero.
    u=tb-((tb-tc)*q-(tb-ta)*r)...
        /(2.0*sign(max(abs(q-r),TINY),q-r));
    ulim=tb+GLIMIT*(tc-tb);

    if ((tb-u)*(u-tc) > 0.0) % u is between tb and tc: try it.
        fu=myfeval(f,v0+u*vec,data);
        if (fu < fc)  % Got a minimum between b and c
            ta=tb; tb=u;
            fa=fb; fb=fu;
            done=1;
            continue;
        elseif (fu > fb) % Got a minimum between a and u
            tc=u;
            fc=fu;
            done=1;
            continue;
        end
        u=tc+GOLD*(tc-tb); % Parabolic fit was of no use. Use default magnification.
        fu=myfeval(f,v0+u*vec,data);
    elseif ((tc-u)*(u-ulim) > 0.0) % Parabolic fit is between c and its allowed limit
        fu=myfeval(f,v0+u*vec,data);
        if (fu < fc)
            tb=tc; tc=u; u=tc+GOLD*(tc-tb);
            fb=fc; fc=fu; fu=myfeval(f,v0+u*vec,data);
        end
    elseif ((u-ulim)*(ulim-tc) >= 0.0) % Limit parabolic u to maximum allowed value
        u=ulim;
        fu=myfeval(f,v0+u*vec,data);
    else
        u=tc+GOLD*(tc-tb); % Reject parabolic u, use default magnification.
        fu=myfeval(f,v0+u*vec,data);
    end
    
    ta=tb; tb=tc; tc=u; % Eliminate oldest point and continue.
    fa=fb; fb=fc; fc=fu;
    done=(fc >= fb);
end

% convert from one-dimensional bracketing parameter along the line from ax
% to bx, to bracketing points in n-dimensions.
ax=v0+ta*vec;
bx=v0+tb*vec;
cx=v0+tc*vec;

if DEBUG
    Nt=100;
    t=linspace(ta,tc,Nt);
    vals=zeros(Nt,1);
    for j=1:Nt
        vals(j)=myfeval(f,v0+t(j)*vec,data);
    end
    plot(t,vals);
    hold on;
    scatter(ta,fa,15,'r','filled');
    scatter(tb,fb,15,'r','filled');
    scatter(tc,fc,15,'r','filled');
    text(ta,fa+0.1,'a');
    text(tb,fb+0.1,'b');
    text(tc,fc+0.1,'c');
    hold off;
end

return;


function  c=sign(a,b)
if b>=0
    c=abs(a);
else
    c=-abs(a);
end
    

%% testing
%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) x*(x-2),[],1)
%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) x*(x+2),[],1)
%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) cos(x),[],1)
%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(10*x),[],1)
%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(10*x)+cos(x),[],1)
%[ax,bx,cx,fa,fb,fc]=mnbrak(0,1,@(x) sin(x),[],1)
%[ax,bx,cx,fa,fb,fc,err]=mnbrak(0,1,@(x,j) besselj(0,x.^j),0.5,1)
%[ax,bx,cx,fa,fb,fc,err]=mnbrak([0 0 0],[1 1 1],@(x) besselj(0,norm(x)),[],1)
%[ax,bx,cx,fa,fb,fc,err]=mnbrak([0 0 0],[1 1 1],@(x) sum((x-[1 1 1]).*x),[],1)
%[ax,bx,cx,fa,fb,fc,err]=mnbrak(0,1,@(x) -(x),[],1) %err=1

