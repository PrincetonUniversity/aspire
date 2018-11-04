function plot_view_dirs(viewingdirs,color)

thetas = atan2d(viewingdirs(2,:),viewingdirs(1,:));
phis = acosd(viewingdirs(3,:));

folded_phis = phis;
idx = find(folded_phis>90);
folded_phis(idx) = 180-folded_phis(idx);
assert(all(folded_phis<=90))

style = sprintf('o%s',color);
polar(thetas/180*pi,folded_phis,style);

end