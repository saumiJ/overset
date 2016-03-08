function isWorking = overset_test_grid()
% basic test function for overset grid
clear
clc
close all

% base-grid params
bg_id = 0;
bg_nx = 26; bg_ny = 26;
bg_dx = 1; bg_dy = 1;

% sub-grid-1 params
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

% sub-grid-2 params
sg2_id = 2;
sg2_nx = 8; sg2_ny = 8;
sg2_dx = 1; sg2_dy = 1;
sg2_center_init = [10 10];
sg2_angle_init = 3*pi/8;

% create base- and sub-grids
oBG = overset_base_grid('T_b', bg_id, bg_nx, bg_ny, bg_dx, bg_dy, {});
oSG1 = overset_sub_grid('T_sub1', sg1_id, sg1_nx, sg1_ny, sg1_dx, sg1_dy, sg1_list_voidBoundaryPointList, sg1_center_init, sg1_angle_init);
oSG2 = overset_sub_grid('T_sub2', sg2_id, sg2_nx, sg2_ny, sg2_dx, sg2_dy, {}, sg2_center_init, sg2_angle_init);

% choose number of grids to work with
n_grids = 2;
if n_grids == 2
    grids = cell(1, 2);
    grids{1} = oBG;
    grids{2} = oSG1;
else if n_grids == 3
    grids = cell(1, 3);
    grids{1} = oBG;
    grids{2} = oSG1;
    grids{3} = oSG2;
    else
        disp 'Check n_grids!'
        exit(1)
    end
end

% create composite grid
overlap = 1;
oCG = overset_composite_grid('T_comp', grids, overlap);

% interpolate values
oCG.interpolate();

% plot grid
fig_grid = figure();
oCG.display_grid(fig_grid);

% print data
fig_data = figure();

% solve poisson problem!
T_vector = [100 300];
isOuterNeumann = true;
n_iter = 15;
oCG = overset_steady_poisson_problem(oCG, T_vector, isOuterNeumann, n_iter);

oCG.display_data(fig_data);

isWorking = 1;
end
