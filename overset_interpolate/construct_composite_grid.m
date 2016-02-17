function [] = construct_composite_grid(grids)
% constructs composite grid from 'grids'
% note: 'grids' must be arranged such that the farthest one will be the
% topmost, i.e. it will not undergo cell removal (unless points lie outside
% domain)

% ref: Cesshire et. al. (1990). Composite Overlapping Meshes for the 
% Solution of Partial Differential Equations

n_grids = size(grids, 2);

% STEP 1: assign flags to n_grids
for k = 1: n_grids
    for i = 1: grids{k}.ny
        for j = 1: grids{k}.nx
            grids{k}.flag(i, j) = n_grids;
        end
    end
end

% STEP 2: eliminate non-boundary points close to boundary points of other grids

% repeatedly needed quantities
boundary_indices_i = cell(n_grids, 1);
boundary_indices_j = cell(n_grids, 1);
boundaryButOne_indices_i = cell(n_grids, 1);
boundaryButOne_indices_j = cell(n_grids, 1);
voidBoundary_indices_i = cell(n_grids, 1);
voidBoundary_indices_j = cell(n_grids, 1);
polygon_x = cell(n_grids, 1);
polygon_y = cell(n_grids, 1);
global_coords = cell(n_grids, 1);
for k = 1: n_grids
    % get boundary indices
    boundary_indices_i{k} = [ones(1, grids{k}.nx-1) ...
        1: grids{k}.ny-1 ...
        grids{k}.ny*ones(1, grids{k}.nx-1) ...
        grids{k}.nx: -1: 2];
    
    boundary_indices_j{k} = [1: grids{k}.nx-1 ...
        grids{k}.nx*ones(1, grids{k}.ny-1) ...
        grids{k}.ny: -1: 2 ...
        ones(1, grids{k}.ny-1)];
    
    voidBoundary_indices_i{k} = [];
    voidBoundary_indices_j{k} = [];
    for i = 1: grids{k}.ny
        for j = 1: grids{k}.nx
            if grids{k}.isVoidBoundary(i, j) == 1
                voidBoundary_indices_i{k} = [voidBoundary_indices_i{k} i];
                voidBoundary_indices_j{k} = [voidBoundary_indices_j{k} j];
            end
        end
    end
    % get boundary-but-one indices
    boundaryButOne_indices_i{k} = [2*ones(1, grids{k}.nx-3) ...
        2: grids{k}.ny-3 ...
        (grids{k}.ny-1)*ones(1, grids{k}.nx-3) ...
        grids{k}.nx-1: -1: 3];
    
    boundaryButOne_indices_j{k} = [2: grids{k}.nx-3 ...
        (grids{k}.nx-1)*ones(1, grids{k}.ny-3) ...
        grids{k}.ny-1: -1: 3 ...
        2*ones(1, grids{k}.ny-3)];
        
    % global coords for grid k
    global_coords{k} = grids{k}.get_global_coords();
    
    % construct polygon points
    polygon_x{k} = [global_coords{k}(1, 1, 2) ...
                global_coords{k}(1, end, 2) ...
                global_coords{k}(end, end, 2) ...
                global_coords{k}(end, 1, 2)];
            
    polygon_y{k} = [global_coords{k}(1, 1, 1) ...
                global_coords{k}(1, end, 1) ...
                global_coords{k}(end, end, 1) ...
                global_coords{k}(end, 1, 1)];
end

for k = 1: n_grids
    for kd = 1: n_grids
        if kd ~= k
                        
            % make a list of void-boundary points in k
            k_voidBoundary_points = zeros(length(voidBoundary_indices_i{k}), 2);
            for l = 1: length(voidBoundary_indices_i{k})
                k_voidBoundary_points(l, 1) = global_coords{k}(voidBoundary_indices_i{k}(l), voidBoundary_indices_j{k}(l), 1);
                k_voidBoundary_points(l, 2) = global_coords{k}(voidBoundary_indices_i{k}(l), voidBoundary_indices_j{k}(l), 2);
            end
            
            % check which boundary points of k lie in kd
            is_point_in_kd = inpolygon(k_voidBoundary_points(:, 2), k_voidBoundary_points(:, 1), polygon_x{kd}, polygon_y{kd});
            pick_points_at = find(is_point_in_kd==1);
            voidBoundary_indices_i_in_kd = zeros(1, length(pick_points_at));
            voidBoundary_indices_j_in_kd = zeros(1, length(pick_points_at));
            for l = 1: length(pick_points_at)
                voidBoundary_indices_i_in_kd(l) = boundary_indices_i{k}(pick_points_at(l));
                voidBoundary_indices_j_in_kd(l) = boundary_indices_j{k}(pick_points_at(l));
            end

            % loop over kd-points, find distances between
            % this and k-voidBoundary pts
            if ~isempty(k_voidBoundary_points)
                for l = 1: length(k_voidBoundary_points)
                    min_dist_indices = []; min_euclidean_dist = []; found = false;
                    for i = 1: grids{kd}.ny
                        for j = 1: grids{kd}.nx
                            kd_point_to_check = [global_coords{kd}(i, j, 1) global_coords{kd}(i, j, 2)];
                            manhattan_dist = abs(kd_point_to_check(1, :) - k_voidBoundary_points(l, :));
                            if (manhattan_dist(1) < grids{kd}.dy/2) && (manhattan_dist(2) < grids{kd}.dx/2)
                                if found == false
                                    min_dist_indices = [i j];
                                    min_euclidean_dist = sqrt(sum(bsxfun(@minus, kd_point_to_check, k_voidBoundary_points(l, :)).^2,2));
                                    found = true;
                                else
                                    euclidean_dist = sqrt(sum(bsxfun(@minus, kd_point_to_check, k_voidBoundary_points(l, :)).^2,2));
                                    if euclidean_dist < min_euclidean_dist
                                        min_dist_indices = [i j];
                                        min_euclidean_dist = euclidean_dist;
                                    end
                                end
                            end
                        end
                    end
                    % if closer to void-boundary, set that point's flag to 0
                    if found
                        grids{kd}.flag(min_dist_indices(1, 1), min_dist_indices(1, 2)) = 0;
                    end
                end
            end
        end
    end
end

% STEP 3: Find and mark points with grids they can be interpolated from. If
% not, mark as exterior or interior
isChanging = true;
% loop while atleast one flag changes
while isChanging
    isChanging = false;
    for k = 1: n_grids
        for i = 1: grids{k}.ny
            for j = 1: grids{k}.nx
                kd = grids{k}.flag(i, j);
                if kd ~= 0 % for all non-exterior grid points
                    isWhileLoopDone = false;
                    while ~isWhileLoopDone
                        isValidPoint = false;
                        if kd == k
                            % Check if grids{k}.grid_coords(i, j, :) is a valid
                            % discretization point.
                            for p = 1: n_grids
                                if grids{p}.num_void_polygons > 0
                                    l = 1; 
                                    while (l <= grids{p}.num_void_polygons) && ~isValidPoint
                                        [poly_x, poly_y] = grids{p}.get_void_polygon(l);
                                        isValidPoint = ~inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), poly_x, poly_y);
                                        l = l + 1;
                                    end
                                end
                            end
                            if ~isValidPoint
                                grids{k}.flag(i, j) = grids{k}.flag(i, j) - 1;
                                isChanging = true;
                            end
                        else
                            % check if point can be interpolated from kd grid
                            % TODO: linear interpolation for now. Extend to
                            % arbitrary order of interpolation
                            iskPointInkd = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_x{kd}, polygon_y{kd});
                            if iskPointInkd  % if point in kd
                                % check if it falls in a void
                                iskPointInkdVoid = false;
                                l = 1;
                                while (l <= grids{kd}.num_void_polygons)
                                    [poly_x, poly_y] = grids{kd}.get_void_polygon(l);
                                    iskPointInkdVoid = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), poly_x, poly_y);
                                    l = l + 1;
                                end
                                if iskPointInkdVoid % if in void, cannot interpolate
                                    grids{k}.flag(i, j) = grids{k}.flag(i, j) - 1;
                                    isChanging = true;
                                end
                            else % if point cannot be interpolated from kd grid
                                grids{k}.flag(i, j) = grids{k}.flag(i, j) - 1;
                                isChanging = true;
                            end
                        end
                        kd = kd - 1;
                        if kd == 0 || isValidPoint
                            isWhileLoopDone = true;
                            isChanging = false;
                        end
                    end
                end
            end
        end
    end
end

end