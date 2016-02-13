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

% STEP 2: eliminate boundary points close to boundary points of other grids
boundary_indices_i = cell(n_grids, 1);
boundary_indices_j = cell(n_grids, 1);
boundaryButOne_indices_i = cell(n_grids, 1);
boundaryButOne_indices_j = cell(n_grids, 1);
global_coords = cell(n_grids, 1);
for k = 1: n_grids
    boundary_indices_i{k} = [ones(1, grids{k}.nx-1) ...
        1: grids{k}.ny-1 ...
        grids{k}.ny*ones(1, grids{k}.nx-1) ...
        grids{k}.nx: -1: 2];
    
    boundary_indices_j{k} = [1: grids{k}.nx-1 ...
        grids{k}.nx*ones(1, grids{k}.ny-1) ...
        grids{k}.ny: -1: 2 ...
        ones(1, grids{k}.ny-1)];
    
    boundaryButOne_indices_i{k} = [2*ones(1, grids{k}.nx-3) ...
        2: grids{k}.ny-3 ...
        (grids{k}.ny-1)*ones(1, grids{k}.nx-3) ...
        grids{k}.nx-1: -1: 3];
    
    boundaryButOne_indices_j{k} = [2: grids{k}.nx-3 ...
        (grids{k}.nx-1)*ones(1, grids{k}.ny-3) ...
        grids{k}.ny-1: -1: 3 ...
        2*ones(1, grids{k}.ny-3)];
    
    global_coords{k} = grids{k}.get_global_coords();
    
end

for k = 1: n_grids
    for kd = 1: n_grids
        if kd ~= k
            % construct polygon points of grid kd
            polygon_kd_x = [global_coords{kd}(1, 1, 2) ...
                global_coords{kd}(1, end, 2) ...
                global_coords{kd}(end, end, 2) ...
                global_coords{kd}(end, 1, 2)];
            polygon_kd_y = [global_coords{kd}(1, 1, 1) ...
                global_coords{kd}(1, end, 1) ...
                global_coords{kd}(end, end, 1) ...
                global_coords{kd}(end, 1, 1)];

            % make a list of boundary points in k
            k_boundary_points = zeros(length(boundary_indices_i{k}), 2);
            for l = 1: length(boundary_indices_i{k})
                k_boundary_points(l, 1) = global_coords{k}(boundary_indices_i{k}(l), boundary_indices_j{k}(l), 1);
                k_boundary_points(l, 2) = global_coords{k}(boundary_indices_i{k}(l), boundary_indices_j{k}(l), 2);
            end
            
            % check which boundary points of k lie in kd
            is_point_in_kd = inpolygon(k_boundary_points(:, 2), k_boundary_points(:, 1), polygon_kd_x, polygon_kd_y);
            pick_points_at = find(is_point_in_kd==1);
            boundary_indices_i_in_kd = zeros(1, length(pick_points_at));
            boundary_indices_j_in_kd = zeros(1, length(pick_points_at));
            
            for l = 1: length(pick_points_at)
                boundary_indices_i_in_kd(l) = boundary_indices_i{k}(pick_points_at(l));
                boundary_indices_j_in_kd(l) = boundary_indices_j{k}(pick_points_at(l));
            end
            
            % check if kd point that is closest to shortlisted k
            % boundary point is a kd boundary point
            
            % make a list of boundaryButOne points in kd
            kd_boundaryButOne_points = zeros(length(boundaryButOne_indices_i{kd}), 2);
            for l = 1: length(boundaryButOne_indices_i{kd})
                kd_boundaryButOne_points(l, 1) = global_coords{kd}(boundaryButOne_indices_i{kd}(l), boundaryButOne_indices_j{kd}(l), 1);
                kd_boundaryButOne_points(l, 2) = global_coords{kd}(boundaryButOne_indices_i{kd}(l), boundaryButOne_indices_j{kd}(l), 2);
            end

            % make a list of boundary points in kd
            kd_boundary_points = zeros(length(boundary_indices_i{kd}), 2);
            for l = 1: length(boundary_indices_i{kd})
                kd_boundary_points(l, 1) = global_coords{kd}(boundary_indices_i{kd}(l), boundary_indices_j{kd}(l), 1);
                kd_boundary_points(l, 2) = global_coords{kd}(boundary_indices_i{kd}(l), boundary_indices_j{kd}(l), 2);
            end

            % loop over k-boundary-points-in-kd, find distances between
            % this and (a) kd-boundaryButOne pts and kd-boundary pts
            for l = 1: length(boundary_indices_i_in_kd)
                k_point_to_check = [global_coords{k}(boundary_indices_i_in_kd(l), boundary_indices_j_in_kd(l), 1) ...
                    global_coords{k}(boundary_indices_i_in_kd(l), boundary_indices_j_in_kd(l), 2)];
                
                distances_boundaryButOne = sqrt(sum(bsxfun(@minus, k_point_to_check, kd_boundaryButOne_points).^2,2));
                distances_boundary       = sqrt(sum(bsxfun(@minus, k_point_to_check, kd_boundary_points).^2,2));
                
                closest_boundary_point_indices = [boundary_indices_i{kd}(distances_boundary==min(distances_boundary)), ...
                    boundary_indices_j{kd}(distances_boundary==min(distances_boundary))];
                
                % if closer to boundary, set that point's flag to 0
                if min(distances_boundary) < min(distances_boundaryButOne)
                    grids{kd}.flag(closest_boundary_point_indices(1), closest_boundary_point_indices(2)) = 0;
                end
            end
            
        end
    end
end

% STEP 3: 

end