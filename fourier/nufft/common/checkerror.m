function checkerror(n,f1,f2,alpha,epsilon,t1,t2)
%
% Check that the vectors f1 and f2 are close enough.
% The vectors are close enough if 
%    max(f1-f2)/sum(abs(alpha))<eps
%
% t1 and t2 (optional) are the times required to compute f1 and f2.
% Typically, f1 should correspond to the slower algorithm.
%
% Yoel Shkolnisky, January 2010.

d=max(f1(:)-f2(:))/sum(abs(alpha(:)));
d=abs(d);

if nargin<6
    % No timings are given
    if abs(d)<epsilon
        str=srpintf('n=%6d  OK    d=%+7.4e\n',n,d);
    else
        str=sprintf('n=%6d  ERROR     d=%+7.4e    norm(f1)=%7.4e    norm(f2)=%7.4e\n',n,d,norm(f1(:)),norm(f2(:)));
    end
else
    % Timings are given
    if t1<1.0e-12
        t1=1;        
    end
    if t2<1.0e-12
        t2=1;
    end
    
    if abs(d)<epsilon
        str=sprintf('n=%6d  OK    d=%+7.4e    t1/t2=%7.4f\n',n,d,t1/t2);
    else
        str=sprintf('n=%6d  ERROR     d=%+7.4e    norm(f1)=%7.4e     norm(f2)=%7.4e\n',n,d,norm(f1(:)),norm(f2(:))); % XXX fix all tests to this form.
    end
end

fprintf(str);


