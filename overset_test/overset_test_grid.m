function isWorking = overset_test_grid()
% basic test function for overset grid
clear
clc
close all

bg_id = 0;
bg_nx = 26; bg_ny = 26;
bg_dx = 1; bg_dy = 1;

sg1_id = 1;
sg1_nx = 15; sg1_ny = 15;
sg1_dx = 1; sg1_dy = 1;
sg1_center_init = [13 13];
sg1_angle_init = pi/6;
sg1_number_of_voids = 1;
sg1_list_voidBoundaryPointList = cell(1, sg1_number_of_voids);
for k = 1: sg1_number_of_voids
    sg1_list_voidBoundaryPointList{k} = [6 6 6 6 6  6  7  8  9  10 11 11 11 11 11 11 10 9 8 7;...
                                         6 7 8 9 10 11 11 11 11 11 11 10 9  8  7  6  6  6 6 6];
end

% sg2_id = 1;
% sg2_nx = 6; sg2_ny = 6;
% sg2_dx = 0.75; sg2_dy = 0.75;
% sg2_center_init = [3 3];
% sg2_angle_init = 2*pi/6;

oBG = overset_base_grid('T_b', bg_id, bg_nx, bg_ny, bg_dx, bg_dy, {});
oSG1 = overset_sub_grid('T_sub1', sg1_id, sg1_nx, sg1_ny, sg1_dx, sg1_dy, sg1_list_voidBoundaryPointList, sg1_center_init, sg1_angle_init);
% oSG2 = overset_sub_grid('T_sub2', sg2_id, sg2_nx, sg2_ny, sg2_dx, sg2_dy, sg2_center_init, sg2_angle_init);

fig_grid = figure();
% oSG2.display_grid(fig_grid);


oBG.display_grid(fig_grid);
oSG1.display_grid(fig_grid);
construct_composite_grid({oBG, oSG1});
% construct_composite_grid({oBG, oSG1, oSG2});

isWorking = 1;
end