function oCG = overset_steady_poisson_problem(oCG, T_vector, isOuterNeumann, n_iter, isPlotSaved, problemName)
% 2D steady-state poisson problem over an overset composite grid

disp '--------------------------------'
disp 'overset transient poisson solver'
disp '--------------------------------'

A = construct_poisson_matrices(oCG, isOuterNeumann);
b = construct_right_hand_sides(oCG, isOuterNeumann, T_vector);

fig_grid = figure();
fig_data = figure();

for iter = 1: n_iter
    disp (strcat('overset: solver | -> iterartion: ', num2str(iter)));
    T_sol = solve(A, b);
    oCG = map_results_to_grid(oCG, T_sol);
    oCG.interpolate();
    b = construct_right_hand_sides(oCG, true, T_vector);
end

oCG.display_data(fig_data);
if isPlotSaved
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 7 7]);
    saveas(figure(fig_data), strcat(problemName, '_data.tiff'));
end
oCG.display_grid(fig_grid);
if isPlotSaved
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 7 7]);
    saveas(figure(fig_grid), strcat(problemName, '_grid.tiff'));
end

end

function A = construct_poisson_matrices(oCG, isOuterNeumann)
% function to construct the solution matrices for the Poisson problem on oCG

A = cell(oCG.n_grids, 1);

% ---------------------
% ---------------------
%      BASE GRID
% ---------------------
% ---------------------

k = 1;

nx = oCG.grids{k}.nx;
ny = oCG.grids{k}.ny;
dx = oCG.grids{k}.dx;
dy = oCG.grids{k}.dy;

i_sp = zeros(1, sqrt(nx*ny)*(ny * nx));
j_sp = zeros(1, sqrt(nx*ny)*(ny * nx));
a_sp = zeros(1, sqrt(nx*ny)*(ny * nx));

ctr = 1;
% ---------------------
% fill corner elements
% ---------------------
i_holder = [
    1 
    nx 
    nx*(ny - 1) + 1 
    ny * nx
    ]; 
j_holder = i_holder;
a_holder = [1 1 1 1];
for l = 1: 4
    i_sp(ctr) = i_holder(l);
    j_sp(ctr) = j_holder(l);
    a_sp(ctr) = a_holder(l);
    ctr = ctr + 1;
end

% ---------------------
% fill lower row
% ---------------------
i_holder = zeros(1, 2*(nx-2));
j_holder = zeros(1, 2*(nx-2));
a_holder = zeros(1, 2*(nx-2));
lim = 0;
if isOuterNeumann
    lim = lim + 2*(nx-2);
    for l = 2: 2: 2*(nx-2)
        i_holder(l-1) = l2Tol1(l);
        j_holder(l-1) = l2Tol1(l);
        a_holder(l-1) = 1;
        
        i_holder(l) = l2Tol1(l);
        j_holder(l) = l2Tol1(l) + nx;
        a_holder(l) = -1;
    end
else
    lim = lim + (nx-2);
    for l = 2: (nx-1)
        i_holder(l-1) = l;
        j_holder(l-1) = l;
        a_holder(l-1) = 1;
    end
end

for l = 1: lim
    i_sp(ctr) = i_holder(l);
    j_sp(ctr) = j_holder(l);
    a_sp(ctr) = a_holder(l);
    ctr = ctr + 1;
end

% ---------------------
% fill upper row
% ---------------------
i_holder = zeros(1, 2*(nx-2));
j_holder = zeros(1, 2*(nx-2));
a_holder = zeros(1, 2*(nx-2));
lim = 0;
if isOuterNeumann
    lim = lim + 2*(nx-2);
    for l = 2: 2: 2*(nx-2)
        i_holder(l-1) = l2Tol1(l) + (ny - 1) * nx;
        j_holder(l-1) = l2Tol1(l) + (ny - 1) * nx;
        a_holder(l-1) = 1;

        i_holder(l) = l2Tol1(l) + (ny - 1) * nx;
        j_holder(l) = l2Tol1(l) + (ny - 2) * nx;
        a_holder(l) = -1;
    end
else
    lim = lim + nx-2;
    for l = 2: nx-1
        i_holder(l-1) = l + (ny - 1) * nx;
        j_holder(l-1) = l + (ny - 1) * nx;
        a_holder(l-1) = 1;
    end
end

for l = 1: lim
    i_sp(ctr) = i_holder(l);
    j_sp(ctr) = j_holder(l);
    a_sp(ctr) = a_holder(l);
    ctr = ctr + 1;
end

% ---------------------
% fill left row
% ---------------------
i_holder = zeros(1, 2*(ny-2));
j_holder = zeros(1, 2*(ny-2));
a_holder = zeros(1, 2*(ny-2));
lim = 0;
if isOuterNeumann
    lim = lim + 2*(ny-2);
    for l = 2: 2: 2*(ny-2)
        i_holder(l-1) = 1 + (l2Tol1(l) -1) * nx;
        j_holder(l-1) = 1 + (l2Tol1(l) -1) * nx;
        a_holder(l-1) = 1;
        
        i_holder(l) = 1 + (l2Tol1(l) -1) * nx;
        j_holder(l) = 1 + (l2Tol1(l) -1) * nx + 1;
        a_holder(l) = -1;
    end
else
    lim = lim + ny-2;
    for l = 2: ny-1
        i_holder(l-1) = 1 + (l -1) * nx;
        j_holder(l-1) = 1 + (l -1) * nx;
        a_holder(l-1) = 1;
    end
end

for l = 1: lim
    i_sp(ctr) = i_holder(l);
    j_sp(ctr) = j_holder(l);
    a_sp(ctr) = a_holder(l);
    ctr = ctr + 1;
end

% ---------------------
% fill right row
% ---------------------
i_holder = zeros(1, 2*(ny-2));
j_holder = zeros(1, 2*(ny-2));
a_holder = zeros(1, 2*(ny-2));
lim = 0;
if isOuterNeumann
    lim = lim + 2*(ny-2);
    for l = 2: 2: 2*(ny-2)
        i_holder(l-1) = l2Tol1(l) * nx;
        j_holder(l-1) = l2Tol1(l) * nx;
        a_holder(l-1) = 1;

        i_holder(l) = l2Tol1(l) * nx;
        j_holder(l) = l2Tol1(l) * nx - 1;
        a_holder(l) = -1;
    end
else
    lim = lim + ny-2;
    for l = 2: ny-1
        i_holder(l-1) = l * nx;
        j_holder(l-1) = l * nx;
        a_holder(l-1) = 1;
    end
end

for l = 1: lim
    i_sp(ctr) = i_holder(l);
    j_sp(ctr) = j_holder(l);
    a_sp(ctr) = a_holder(l);
    ctr = ctr + 1;
end

% ---------------------
% fill inner cells
% ---------------------
for i = 2: ny-1
    for j = 2: nx-1
        switch oCG.grids{k}.flag(i, j)
            case 0
                i_sp(ctr) = (i-1) * nx + j;
                j_sp(ctr) = (i-1) * nx + j;
                a_sp(ctr) = 1;
                ctr = ctr + 1;
            case k
                i_holder = [
                    (i-1) * nx + j
                    (i-1) * nx + j
                    (i-1) * nx + j
                    (i-1) * nx + j
                    (i-1) * nx + j
                    ]; 
                j_holder = [
                    (i-1) * nx + j
                    (i-1) * nx + j - 1;
                    (i-1) * nx + j + 1;
                    (i-2) * nx + j;
                    (i-0) * nx + j;
                    ];
                a_holder = [
                    2 * (dx^2 + dy^2)
                    -dx^2
                    -dx^2
                    -dy^2
                    -dy^2
                    ];
                for l = 1: length(a_holder)
                        i_sp(ctr) = i_holder(l);
                        j_sp(ctr) = j_holder(l);
                        a_sp(ctr) = a_holder(l);
                        ctr = ctr + 1;
                end
            otherwise
                i_sp(ctr) = (i-1) * nx + j;
                j_sp(ctr) = (i-1) * nx + j;
                a_sp(ctr) = 1;
                ctr = ctr + 1;
        end
    end
end

A{k} = sparse(i_sp(1: ctr-1), j_sp(1: ctr-1), a_sp(1: ctr-1));

% ---------------------
% loop over all grids
% ---------------------
for k = 2: oCG.n_grids
    nx = oCG.grids{k}.nx;
    ny = oCG.grids{k}.ny;
    dx = oCG.grids{k}.dx;
    dy = oCG.grids{k}.dy;

    i_sp = zeros(1, sqrt(nx*ny)*(ny * nx));
    j_sp = zeros(1, sqrt(nx*ny)*(ny * nx));
    a_sp = zeros(1, sqrt(nx*ny)*(ny * nx));
    
    ctr = 1;
    
    % ---------------------
    % fill inner cells
    % ---------------------
    for i = 1: ny
        for j = 1: nx
            switch oCG.grids{k}.flag(i, j)
                case 0
                    i_sp(ctr) = (i-1) * nx + j;
                    j_sp(ctr) = (i-1) * nx + j;
                    a_sp(ctr) = 1;
                    ctr = ctr + 1;
                case k
                    if oCG.grids{k}.isVoidBoundary(i, j)
                        i_sp(ctr) = (i-1) * nx + j;
                        j_sp(ctr) = (i-1) * nx + j;
                        a_sp(ctr) = 1;
                        ctr = ctr + 1;
                    else
                        i_holder = [
                            (i-1) * nx + j
                            (i-1) * nx + j
                            (i-1) * nx + j
                            (i-1) * nx + j
                            (i-1) * nx + j
                            ]; 
                        j_holder = [
                            (i-1) * nx + j
                            (i-1) * nx + j - 1;
                            (i-1) * nx + j + 1;
                            (i-2) * nx + j;
                            (i-0) * nx + j;
                            ];
                        a_holder = [
                            2 * (dx^2 + dy^2)
                            -dx^2
                            -dx^2
                            -dy^2
                            -dy^2
                            ];
                        for l = 1: length(a_holder)
                                i_sp(ctr) = i_holder(l);
                                j_sp(ctr) = j_holder(l);
                                a_sp(ctr) = a_holder(l);
                                ctr = ctr + 1;
                        end
                    end
                otherwise
                    i_sp(ctr) = (i-1) * nx + j;
                    j_sp(ctr) = (i-1) * nx + j;
                    a_sp(ctr) = 1;
                    ctr = ctr + 1;
            end
        end
    end
    A{k} = sparse(i_sp(1: ctr-1), j_sp(1: ctr-1), a_sp(1: ctr-1));
end

end

function b = construct_right_hand_sides(oCG, isOuterNeumann, T_vector)
% function to compute RHS vectors

b = cell(oCG.n_grids, 1);

% ---------------------
% ---------------------
%      BASE GRID
% ---------------------
% ---------------------

k = 1;

nx = oCG.grids{k}.nx;
ny = oCG.grids{k}.ny;

b{k} = zeros((ny * nx), 1);


% ---------------------
% fill corner elements
% ---------------------
i_holder = [
    1 
    nx 
    nx*(ny - 1) + 1 
    ny * nx
    ]; 
a_holder = [0 0 0 0];
for l = 1: 4
    b{k}(i_holder(l)) = a_holder(l);
end

% ---------------------
% fill lower row
% ---------------------
i_holder = zeros(1, 1*(nx-2));
a_holder = zeros(1, 1*(nx-2));
for l = 1: (nx-2)
    i_holder(l) = l;
    if isOuterNeumann
        a_holder(l) = 0;
    else
        a_holder(l) = T_vector(1);
    end
end
for l = 1: length(a_holder)
    b{k}(i_holder(l)) = a_holder(l);
end

% ---------------------
% fill upper row
% ---------------------
i_holder = zeros(1, 1*(nx-2));
a_holder = zeros(1, 1*(nx-2));
for l = 1: (nx-2)
    i_holder(l) = nx * (ny-1) + 1 + l;
    if isOuterNeumann
        a_holder(l) = 0;
    else
        a_holder(l) = T_vector(1);
    end
end
for l = 1: length(a_holder)
    b{k}(i_holder(l)) = a_holder(l);
end

% ---------------------
% fill left row
% ---------------------
i_holder = zeros(1, 1*(ny-2));
a_holder = zeros(1, 1*(ny-2));
for l = 1: (ny-2)
    i_holder(l) = l*nx + 1;
    if isOuterNeumann
        a_holder(l) = 0;
    else
        a_holder(l) = T_vector(1);
    end
end
for l = 1: length(a_holder)
    b{k}(i_holder(l)) = a_holder(l);
end

% ---------------------
% fill right row
% ---------------------
i_holder = zeros(1, 1*(ny-2));
a_holder = zeros(1, 1*(ny-2));
for l = 1: (ny-2)
    i_holder(l) = (l+1)*nx;
    if isOuterNeumann
        a_holder(l) = 0;
    else
        a_holder(l) = T_vector(1);
    end
end
for l = 1: length(a_holder)
    b{k}(i_holder(l)) = a_holder(l);
end

% ---------------------
% fill inner cells
% ---------------------
for i = 2: ny-1
    for j = 2: nx-1
        ind = (i-1) * nx + j;
        switch oCG.grids{k}.flag(i, j)
            case 0
                b{k}(ind) = 0;
            case k
                b{k}(ind) = 0;
            otherwise
                b{k}(ind) = oCG.grids{k}.val(i, j); % interpolated val
        end
    end
end

% ---------------------
% loop over all grids
% ---------------------
for k = 2: oCG.n_grids
    nx = oCG.grids{k}.nx;
    ny = oCG.grids{k}.ny;
    
    b{k} = zeros((ny * nx), 1);
    
    % ---------------------
    % fill inner cells
    % ---------------------
    for i = 1: ny
        for j = 1: nx
            ind = (i-1) * nx + j;
            switch oCG.grids{k}.flag(i, j)
                case 0
                    b{k}(ind) = oCG.grids{k}.val(i, j); % default void value
                case k
                    if oCG.grids{k}.isVoidBoundary(i, j)
                        b{k}(ind) = T_vector(2); % inner temperature
                    else
                        b{k}(ind) = 0;
                    end
                otherwise
                    b{k}(ind) = oCG.grids{k}.val(i, j); % interpolated val
            end
        end
    end
end
end

function oCG = map_results_to_grid(oCG, T_sol)
% maps result back to composite grid

for k = 1: oCG.n_grids
    for i = 1: oCG.grids{k}.ny
        for j = 1: oCG.grids{k}.nx
            oCG.grids{k}.val(i, j) = T_sol{k}((i-1) * oCG.grids{k}.nx + j);
        end
    end
end

end

function T_sol = solve(A, b)

disp 'overset: solving..'

% function to solve on all grids

n_grids = length(A);
T_sol = cell(n_grids);
for k = 1: n_grids
   T_sol{k} = A{k}\b{k};
end

end

function result = l2Tol1(l)
result = floor(l/2) + mod(l+1,2);
end