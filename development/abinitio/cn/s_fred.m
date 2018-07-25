cryo_abinitio_C4('/home/yoel/scratch/fred/aspire89/averages_nn100_group1.mrc','fred_c4_89_nn100_grp1_2000_linspace.mrc','fred_c4_89_nn100_grp1_2000_linspace.mat',2000);
clear; close all;
cryo_abinitio_C4('/home/yoel/scratch/fred/aspire89/averages_nn100_group2.mrc','fred_c4_89_nn100_grp2_2000_linspace.mrc','fred_c4_89_nn100_grp2_2000_linspace.mat',2000);
clear; close all;

vol1 = ReadMRC('fred_c4_89_nn100_grp1_2000_linspace.mrc');
vol2 = ReadMRC('fred_c4_89_nn100_grp2_2000_linspace.mrc');
pixA = 3.3;
[Rest,estdx,vol2aligned] = cryo_align_densities_C4(vol1,vol2,pixA,1,[],0,50);
[resA,h] = plotFSC(vol1,vol2aligned,0.143,pixA);

% clear
% 
% load fred_c4_89_nn100_grp1_2000_linspace.mat
% cryo_plot_viewing_directions(rots)
