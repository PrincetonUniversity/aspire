function viz_cls(cijs_found,cjis_found,cijs_gt,cjis_gt,n_theta)

cijs_found = cijs_found/n_theta*2*pi;
cjis_found = cjis_found/n_theta*2*pi;
cijs_gt    = cijs_gt/n_theta*2*pi;
cjis_gt    = cjis_gt/n_theta*2*pi;

r_txt_gt    = 1.1;
r_txt_found = 1.2;

figure(10);

subplot(1,2,1);

theta = 0:pi/50:2*pi;
plot(cos(theta),sin(theta));
axis_equal;

for s=1:5
    pos_found = [cos(cijs_found(s)),sin(cijs_found(s))];
    line([0,pos_found(1)],[0,pos_found(2)],'LineWidth',4,'Color','red');
    text(r_txt_found*pos_found(1),r_txt_found*pos_found(2),num2str(s),'Color','red');
    
    pos_gt = [cos(cijs_gt(s)),sin(cijs_gt(s))];
    line([0,pos_gt(1)],[0,pos_gt(2)],'LineWidth',4,'Color','green');
    text(r_txt_gt*pos_gt(1),r_txt_gt*pos_gt(2),num2str(s),'Color','green');
end
axis off;

subplot(1,2,2);

theta = 0:pi/50:2*pi;
plot(cos(theta),sin(theta));
axis_equal;

for s=1:5
    pos_found = [cos(cjis_found(s)),sin(cjis_found(s))];
    line([0,pos_found(1)],[0,pos_found(2)],'LineWidth',4,'Color','red');
    text(r_txt_found*pos_found(1),r_txt_found*pos_found(2),num2str(s),'Color','red');
    
    pos_gt = [cos(cjis_gt(s)),sin(cjis_gt(s))];
    line([0,pos_gt(1)],[0,pos_gt(2)],'LineWidth',4,'Color','green');
    text(r_txt_gt*pos_gt(1),r_txt_gt*pos_gt(2),num2str(s),'Color','green');
end


axis off;