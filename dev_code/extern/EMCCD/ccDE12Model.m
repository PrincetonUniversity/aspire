function y=ccDE12Model(f,theta,p)
% function y=ccDE12Model(f,theta,p)
% Model function for the DE-12 camera, including angle dependence.
fc1=p(1); fc2=p(2); fc3=p(3); a=p(4); a0=p(5); fc0=p(6);
fa0=p(7);
ft0=p(8);
ft0=max(.001,ft0);

ctheta=cos(2*theta);

%     fc0=fc0*(1+ctheta*ft0);
%     fc1=fc1*(1+ctheta*ft1);
%     fc2=fc2*(1+ctheta*ft2);
%     fc3=fc3*(1+ctheta*ft3);
%     a0=a0*(1+ctheta*fa0);
y=a./(1+(f./fc1).^2+(f./fc2).^4+(f./fc3).^6)+a0./(1+(f./fc0).^2);
y=y.*(1+ctheta.*fa0.*(1-1./(1+(f./ft0).^2)));
end
