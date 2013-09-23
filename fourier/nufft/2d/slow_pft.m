function y=slow_pft(x,n_r,n_theta)

if (ndims(x)~=2) || (size(x,1)~=size(x,2))
    error('x must be a 2D square array');
end

dr=2*pi/(2*n_r+1);
dtheta=2*pi/n_theta;
omega=zeros(n_r*n_theta,2);

for k=0:n_theta-1
    for j=0:n_r-1
        omega(k*n_r+j+1,:)=[j*dr*sin(k*dtheta) j*dr*cos(k*dtheta)];
    end
end

y=slow_nufft_t_2d(x,omega);
y=reshape(y,n_r,n_theta);

