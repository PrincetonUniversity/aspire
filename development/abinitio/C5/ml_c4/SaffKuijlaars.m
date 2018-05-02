
function [mesh]=SaffKuijlaars(N)
    
    k=1:N;
    h=-1+2*(k-1)/(N-1);
    theta=acos(h);
    phi=zeros(1,N);
    for k=2:N-1
        phi(k)=mod(phi(k-1)+3.6/(sqrt(N*(1-h(k)^2))),2*pi);
    end
    mesh=(cat(1,sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta)))';
end