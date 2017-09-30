clear;
n_c = 360*180;
n_ic = 600*600*360;

n_tests = 100;
times1 = zeros(1,n_tests);
times2 = zeros(1,n_tests);
for i=1:n_tests
    i
    c = randn(n_c,1);
    ic = ceil(n_c*rand(n_ic,1));
    
    t1 = clock;
    a = c(ic);
    t2 = clock;
    times1(i) = etime(t2,t1);
    
       
    t1 = clock;
    a = unique2all(c,ic);
    t2 = clock;
    times2(i) = etime(t2,t1);
end

mean(times1)
mean(times2)

mean(times1)/mean(times2)

% any(a(:)-aa(:))