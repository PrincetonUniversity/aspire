
function fg=filterGrid(grid,eq_filter_angle)

[eq_idx,~,~]=markEquators(squeeze(grid(:,3,:)),eq_filter_angle);
fg=grid(:,:,~eq_idx);