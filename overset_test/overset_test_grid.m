function isWorking = overset_test_grid()
% basic test function for overset grid
clear
clc
close all

bg_nx = 20; bg_ny = 20;
bg_dx = 1; bg_dy = 1;

sg1_id = 1;
sg1_nx = 20; sg1_ny = 20;
sg1_dx = 0.5; sg1_dy = 0.5;
sg1_center_init = [8.2 8.2];
sg1_angle_init = pi/6;

oBG = overset_base_grid('T_b', bg_nx, bg_ny, bg_dx, bg_dy);
oSG1 = overset_sub_grid('T_sub1', sg1_id, sg1_nx, sg1_ny, sg1_dx, sg1_dy, sg1_center_init, sg1_angle_init);

fig_grid = figure();

oBG.display_grid(fig_grid);
oSG1.display_grid(fig_grid);

isWorking = 1;
end