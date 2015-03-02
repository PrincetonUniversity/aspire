% Test the function fastrotate_ref by comparing it against fastrotate_yaro.
%
% Yoel Shkolnisky, January 2011.

for n=16:2:256
    im=randn(n);
    angle=20;
    % Yaroslavsky's code rotates CW, so take the minus of the angle to get
    % CCW rotation.
    tic; rim1=fastrotate_yaro(im,-angle); t1=toc; 
    tic; rim2=fastrotate_ref(im,angle);  t2=toc;
    
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
    