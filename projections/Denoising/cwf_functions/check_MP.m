function [p,fre]= check_MP(evals, fr, noise_variance, n, nbin)

% Fit the Marcenko Pastur distribution to the given eigenvalues.
% evals: Eigenvalues
% fr: Frequency
% n: Number of data, the n is p/n
% nbin: Number of bins for the histogram
% Tejal Bhamre: Oct 2015

evals=sort(evals,'ascend');
p=size(evals,1);

if fr==0
    y=p/n;
else y=p/(2*n);
end


a = noise_variance*(1-sqrt(y))^2;
b = noise_variance*(1+sqrt(y))^2;
f_MP = @(t) sqrt(max(b-t, 0).*max(t-a, 0) )./(2*pi*y*t*noise_variance);


[nout, xout] = hist(evals, nbin);
hx = xout(2) - xout(1); %step size, used to compute frequency below
extra=(evals(end)-evals(1))/10;
%x1=evals(end)-1;
%x2=evals(1)+1; %two end points
x1=evals(1)-extra;
x2=evals(end)+extra;
xx = x1+hx/2: hx: x2;
fre = f_MP(xx)*hx;


figure,
%axis([0 0.5 0 0.3]);

h = bar(xout, nout/p);
set(h, 'BarWidth', 1, 'FaceColor', 'w', 'EdgeColor', 'b');
hold on;
plot(xx, fre, '--r');
end
