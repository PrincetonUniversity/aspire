% Test the function fastrotate by compring it against fastrotate_ref
%
% Yoel Shkolnisky, January 2011

for n=16:256
    im=randn(n);
    angle=20;      
    t=clock; rim1=fastrotate_ref(im,angle); t1=etime(clock,t);
    t=clock; rim2=fastrotate(im,angle);  t2=etime(clock,t);
    
    diff=rim1-rim2;
    err=norm(diff(:))/norm(im(:));
    
    fprintf('n=%d   speedup=%4.2f  err=%e  ',n,t1/t2,err);
    if err<1.0e-14
        fprintf('OK\n');
    else
        fprintf('ERR\n');
        pause;
    end
end
