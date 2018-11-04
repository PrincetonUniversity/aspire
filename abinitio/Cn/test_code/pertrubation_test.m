function pertrubation_test

ndirs = 10;
refq    = qrand(ndirs);
viewdirs = zeros(3,ndirs);
for i=1:ndirs
    rot = q_to_rot(refq(:,i)).';
    viewdirs(:,i) = rot(:,3);
end

viewdirs_pert = zeros(3,ndirs);
for i=1:ndirs
    view_i = viewdirs(:,i);
    noise = randn(3,1); 
    noise = tand(10)*noise./norm(noise);
    view_pert = view_i + noise;
    viewdirs_pert(:,i) = view_pert./norm(view_pert);
end

figure;
plot_view_dirs(viewdirs,'blue');
hold on;
plot_view_dirs(viewdirs_pert,'red');

end