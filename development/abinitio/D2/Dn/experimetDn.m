
g_y1=axisang2rot([cos(2*pi/3);sin(2*pi/3);0],pi);
g_y2=axisang2rot([cos(4*pi/3);sin(4*pi/3);0],pi);
g_y3=axisang2rot([cos(6*pi/5);sin(6*pi/5);0],pi);
g_y4=axisang2rot([cos(8*pi/5);sin(8*pi/5);0],pi);

g_y1=axisang2rot([cos(2*pi/5);sin(2*pi/5);0],pi);
g_y2=axisang2rot([cos(4*pi/5);sin(4*pi/5);0],pi);
g_y3=axisang2rot([cos(6*pi/5);sin(6*pi/5);0],pi);
g_y4=axisang2rot([cos(8*pi/5);sin(8*pi/5);0],pi);


g_z1=axisang2rot([0,0,1],2*pi/5);
g_z2=axisang2rot([0,0,1],4*pi/5);
g_z3=axisang2rot([0,0,1],6*pi/5);
g_z4=axisang2rot([0,0,1],8*pi/5);


[pf,~]=cryo_pft(projs,size(projs,1),360);
[C]=cryo_calc_proj_corrs(pf(:,:,[19,41]),0,0,1);
[cls,corrs]=cryo_plot_cls_Dn(q_to_rot(q(:,61)),q_to_rot(q(:,62)),3,C);

