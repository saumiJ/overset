classdef overset_composite_grid < handle
% class for composite overset grids
    properties

        % general properties
        name;   % data name

        % data properties
        grids % array of grids
        n_grids; % number of grids
        cell_after_overlap; % first cell after overlap
    
    end
    
    methods
        % constructor
        function obj = overset_composite_grid(name_, grids_, overlap_)
            if length(grids_) < 2
                disp 'ERROR: At least two grids needed to construct a composite grid!'
                exit(1)
            end
            
            gridNames = '';
            for k = 1: length(grids_)
                gridNames = strcat(gridNames, grids_{k}.name, ',');
            end
            disp(strcat('overset: constructing composite grid ', name_, ' from grids: ', gridNames));
            
            obj.name = name_;
            obj.grids = grids_;
            obj.n_grids = length(grids_);
            obj.cell_after_overlap = overlap_ + 1;
            
            obj.construct_composite_grid();
        end
        
        % create grid associations and cuts
        function [] = construct_composite_grid(obj)
            % constructs composite grid from 'grids'
            % note: 'grids' must be arranged such that the farthest one will be the
            % topmost, i.e. it will not undergo cell removal (unless points lie outside
            % domain)

            % ref: Cesshire et. al. (1990). Composite Overlapping Meshes for the 
            % Solution of Partial Differential Equations

            % TODO:
            % - add options for computing cell padding from dx, dy
            % - optimize performance
            % - avoid copying global_coords

            obj.n_grids = size(obj.grids, 2);

            % STEP 1: assign flags to n_grids
            for k = 1: obj.n_grids
                for i = 1: obj.grids{k}.ny
                    for j = 1: obj.grids{k}.nx
                        obj.grids{k}.flag(i, j) = obj.n_grids;
                    end
                end
            end

            % STEP 2: eliminate non-boundary points close to boundary points of other grids

            % repeatedly needed quantities
            polygon_x = cell(obj.n_grids, 1);
            polygon_y = cell(obj.n_grids, 1);
            polygon_but1_x = cell(obj.n_grids, 1);
            polygon_but1_y = cell(obj.n_grids, 1);
            polygon_but2_x = cell(obj.n_grids, 1);
            polygon_but2_y = cell(obj.n_grids, 1);
            global_coords = cell(obj.n_grids, 1);
            voidBoundary_points = cell(obj.n_grids, 1);
            for k = 1: obj.n_grids
                % get void-boundary points
                voidBoundary_points{k} = [];
                for i = 1: obj.grids{k}.num_void_polygons
                    for j = 1: length(obj.grids{k}.void_polygons{i})
                        void_i = obj.grids{k}.void_polygons{i}(1, j);
                        void_j = obj.grids{k}.void_polygons{i}(2, j);
                        point = obj.grids{k}.get_global_coords_at(void_i, void_j);
                        voidBoundary_points{k} = [voidBoundary_points{k}; point];
                    end
                end

                % global coords for grid k
                global_coords{k} = obj.grids{k}.get_global_coords();

                % construct polygon points
                polygon_x{k} = [global_coords{k}(1, 1, 2) ...
                            global_coords{k}(1, end, 2) ...
                            global_coords{k}(end, end, 2) ...
                            global_coords{k}(end, 1, 2)];

                polygon_y{k} = [global_coords{k}(1, 1, 1) ...
                            global_coords{k}(1, end, 1) ...
                            global_coords{k}(end, end, 1) ...
                            global_coords{k}(end, 1, 1)];
                        
                polygon_but1_x{k} = [global_coords{k}(obj.cell_after_overlap,       obj.cell_after_overlap, 2) ...
                                     global_coords{k}(obj.cell_after_overlap,       end-obj.cell_after_overlap+1, 2) ...
                                     global_coords{k}(end-obj.cell_after_overlap+1, end-obj.cell_after_overlap+1, 2) ...
                                     global_coords{k}(end-obj.cell_after_overlap+1, obj.cell_after_overlap, 2)];
                        
                polygon_but1_y{k} = [global_coords{k}(obj.cell_after_overlap,       obj.cell_after_overlap, 1) ...
                                     global_coords{k}(obj.cell_after_overlap,       end-obj.cell_after_overlap+1, 1) ...
                                     global_coords{k}(end-obj.cell_after_overlap+1, end-obj.cell_after_overlap+1, 1) ...
                                     global_coords{k}(end-obj.cell_after_overlap+1, obj.cell_after_overlap, 1)];
                        
                polygon_but2_x{k} = [global_coords{k}(obj.cell_after_overlap+1,     obj.cell_after_overlap+1, 2) ...
                                     global_coords{k}(obj.cell_after_overlap+1,     end-obj.cell_after_overlap, 2) ...
                                     global_coords{k}(end-obj.cell_after_overlap,   end-obj.cell_after_overlap, 2) ...
                                     global_coords{k}(end-obj.cell_after_overlap,   obj.cell_after_overlap+1, 2)];
                        
                polygon_but2_y{k} = [global_coords{k}(obj.cell_after_overlap+1,     obj.cell_after_overlap+1, 1) ...
                                     global_coords{k}(obj.cell_after_overlap+1,     end-obj.cell_after_overlap, 1) ...
                                     global_coords{k}(end-obj.cell_after_overlap,   end-obj.cell_after_overlap, 1) ...
                                     global_coords{k}(end-obj.cell_after_overlap,   obj.cell_after_overlap+1, 1)];
            end

            for k = 1: obj.n_grids
                for kd = 1: obj.n_grids
                    if kd ~= k
                        % loop over kd-points, find distances between
                        % this and k-voidBoundary pts
                        if ~isempty(voidBoundary_points{k})
                            for l = 1: length(voidBoundary_points{k})
                                min_dist_indices = []; min_euclidean_dist = []; found = false;
                                for i = 1: obj.grids{kd}.ny
                                    for j = 1: obj.grids{kd}.nx
                                        kd_point_to_check = [global_coords{kd}(i, j, 1) global_coords{kd}(i, j, 2)];
                                        manhattan_dist = abs(kd_point_to_check(1, :) - voidBoundary_points{k}(l, :));
                                        if (manhattan_dist(1) < obj.grids{kd}.dy/2) && (manhattan_dist(2) < obj.grids{kd}.dx/2)
                                            if ~found
                                                min_dist_indices = [i j];
                                                min_euclidean_dist = sqrt(sum(bsxfun(@minus, kd_point_to_check, voidBoundary_points{k}(l, :)).^2,2));
                                                found = true;
                                            else
                                                euclidean_dist = sqrt(sum(bsxfun(@minus, kd_point_to_check, voidBoundary_points{k}(l, :)).^2,2));
                                                if euclidean_dist < min_euclidean_dist
                                                    min_dist_indices = [i j];
                                                    min_euclidean_dist = euclidean_dist;
                                                end
                                            end
                                        end
                                    end
                                end
                                % if closer to void-boundary, set that kd-point's flag to 0
                                if found
                                    obj.grids{kd}.flag(min_dist_indices(1, 1), min_dist_indices(1, 2)) = 0;
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
                touched = false;
                % loop over grids
                for k = 1: obj.n_grids
                    % loop over all points in the grid
                    for i = 1: obj.grids{k}.ny
                        for j = 1: obj.grids{k}.nx
                            kd = obj.grids{k}.flag(i, j);
                            if kd ~= 0 % for all non-exterior grid points
                                isWhileLoopDone = false;
                                while ~isWhileLoopDone
                                    isValidPoint = false;
                                    if kd == k
                                        % Check if grids{k}.grid_coords(i, j, :) is a valid
                                        % discretization point.
                                        for p = 1: obj.n_grids
                                            if obj.grids{p}.num_void_polygons > 0
                                                l = 1; 
                                                while (l <= obj.grids{p}.num_void_polygons) && ~isValidPoint
                                                    [poly_x, poly_y] = obj.grids{p}.get_void_polygon(l);
                                                    [in, on] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), poly_x, poly_y); % invalid if inside a void_polygon
                                                    isInVoid = in && ~on;
                                                    isValidPoint = ~isInVoid;
                                                    l = l + 1;
                                                end
                                            end
                                        end

                                        % check if on the interface
                                        if (k ~= 1 && isValidPoint)
                                            [in, on] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_x{1}, polygon_y{1});
                                            isInDomain = in && ~on;
                                            [in, on] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_x{k}, polygon_y{k});
                                            isOnInterface = isInDomain && in && on;
                                            isValidPoint = ~isOnInterface;
                                        end

                                        if ~isValidPoint
                                            obj.grids{k}.flag(i, j) = obj.grids{k}.flag(i, j) - 1;
                                            touched = true;
                                        end
                                    else
                                        % check if point can be interpolated from kd grid
                                        % TODO: linear interpolation for now. Extend to
                                        % arbitrary order of interpolation
                                        [in, on] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_x{kd}, polygon_y{kd});
                                        iskPointInkd = in || on;
                                        if iskPointInkd  % if point in kd
                                            % check if it falls in a void
                                            isInVoid = false;
                                            for p = 1: obj.n_grids
                                                if obj.grids{p}.num_void_polygons > 0
                                                    l = 1; 
                                                    while (l <= obj.grids{p}.num_void_polygons) && ~isInVoid
                                                        [poly_x, poly_y] = obj.grids{p}.get_void_polygon(l);
                                                        [in, on] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), poly_x, poly_y); % invalid if inside a void_polygon
                                                        isInVoid = in && ~on;
                                                        l = l + 1;
                                                    end
                                                end
                                            end
                                            isValidPoint = ~isInVoid;
                                            
                                            % check if it needs a kd_boundary point for interpolation
                                            [in, on] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_but1_x{kd}, polygon_but1_y{kd});
                                            isInBut1kd = in || on;

                                            if ~isValidPoint || ~isInBut1kd % if in void, or if needs kd boundary point, cannot interpolate
                                                obj.grids{k}.flag(i, j) = obj.grids{k}.flag(i, j) - 1;
                                                touched = true;
                                            end
                                        else % if point cannot be interpolated from kd grid
                                            obj.grids{k}.flag(i, j) = obj.grids{k}.flag(i, j) - 1;
                                            touched = true;
                                        end
                                    end
                                    kd = kd - 1;
                                    if kd == 0 || isValidPoint
                                        isWhileLoopDone = true;
                                        if ~touched
                                            isChanging = false;
                                        else
                                            isChanging = true;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % STEP 4: Find points on lower grids DEFINITELY needed for interpolation by
            % higher grids

            for k = 1: obj.n_grids
                k_interp_point_count = 0; % counts how many points on k-grid are interpolated from other grids
                k_interp_points = {};
                k_interp_source_ids = {};
                for i = 1: obj.grids{k}.ny
                    for j = 1: obj.grids{k}.nx
                        kd = obj.grids{k}.flag(i, j);
                        if kd < k && kd > 0 % if i, j, k interpolates from kd grid
                            % assuming linear interpolation
                            % TODO: provide arbitrary interpolation order

                            % find the nearest four kd points to k point
                            k_point = [global_coords{k}(i, j, 1) global_coords{k}(i, j, 2)];
                            closest = [inf inf inf inf inf];
                            closest_id = [inf inf inf inf inf];
                            closest_jd = [inf inf inf inf inf];
                            for id = 1: obj.grids{kd}.ny
                                for jd = 1: obj.grids{kd}.nx
                                    %if id == 12 && jd == 8 && i == 8 && j == 7 && k == 3
                                    %    disp 'DEBUG!'
                                    %end
                                    euclidean_dist = sqrt(sum(bsxfun(@minus, k_point, [global_coords{kd}(id, jd, 1) global_coords{kd}(id, jd, 2)]).^2,2));
                                    closest(5) = euclidean_dist; 
                                    closest_id(5) = id;
                                    closest_jd(5) = jd;

                                    [closest, permutor] = sort(closest);
                                    closest_id = closest_id(permutor);
                                    closest_jd = closest_jd(permutor);
                                end
                            end

                            %closest = closest(1: 4);
                            closest_id = closest_id(1: 4);
                            closest_jd = closest_jd(1: 4);

                            % negate kd flag at those 4 pts AND store those 4
                            % pts in correspondence to i j k point for later
                            % use
                            k_interp_point_count = k_interp_point_count + 1;
                            k_interp_points{k_interp_point_count} = [i j];
                            source_ids = zeros(4, 2);
                            for l = 1: 4
                                obj.grids{kd}.flag(closest_id(l), closest_jd(l)) = -abs(obj.grids{kd}.flag(closest_id(l), closest_jd(l)));
                                source_ids(l, :) = [closest_id(l) closest_jd(l)];
                            end
                            k_interp_source_ids{k_interp_point_count} = source_ids;
                        end
                    end
                end
                obj.grids{k}.interp_point_count = k_interp_point_count;
                obj.grids{k}.interp_points = k_interp_points;
                obj.grids{k}.interp_source_ids = k_interp_source_ids;
            end

            % STEP 5: Remove unnecessary interpolation points, change them to
            % discretization points if possible, mark points on higher grids needed for
            % interpolation by lower grids

            for k = 1: obj.n_grids
                k_interp_point_count = 0; % counts how many points on k-grid are interpolated from other grids
                k_interp_points = {};
                k_interp_source_ids = {};
                for i = 1: obj.grids{k}.ny
                    for j = 1: obj.grids{k}.nx
                        % definition of "not needed": this point could be used by
                        % another grid for interpolation, but doesn't need to be
                        kd = abs(obj.grids{k}.flag(i, j));
                        if kd > k
                            % check if point needs interpolation from higher grids
                            [inbB1, onbB1] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_but1_x{kd}, polygon_but1_y{kd});
                            [inbB2, onbB2] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_but2_x{kd}, polygon_but2_y{kd});
                            %[inbB3, onbB3] = inpolygon(global_coords{k}(i, j, 2), global_coords{k}(i, j, 1), polygon_but3_x{kd}, polygon_but3_y{kd});
                            isInBButOne_kd = inbB1 || onbB1;
                            isInBButTwo_kd = inbB2 && ~onbB2;
                            %isInBButThree_kd = inbB3 && ~onbB3;

                            %if isInBButOne_kd && ~isInBButThree_kd
                            if isInBButOne_kd && ~isInBButTwo_kd
                                % find the nearest four kd points to k point
                                k_point = [global_coords{k}(i, j, 1) global_coords{k}(i, j, 2)];
                                closest = [inf inf inf inf inf];
                                closest_id = [inf inf inf inf inf];
                                closest_jd = [inf inf inf inf inf];
                                for id = 1: obj.grids{kd}.ny
                                    for jd = 1: obj.grids{kd}.nx
                                        euclidean_dist = sqrt(sum(bsxfun(@minus, k_point, [global_coords{kd}(id, jd, 1) global_coords{kd}(id, jd, 2)]).^2,2));
                                        closest(5) = euclidean_dist; 
                                        closest_id(5) = id;
                                        closest_jd(5) = jd;

                                        [closest, permutor] = sort(closest);
                                        closest_id = closest_id(permutor);
                                        closest_jd = closest_jd(permutor);
                                    end
                                end

                                %closest = closest(1: 4);
                                closest_id = closest_id(1: 4);
                                closest_jd = closest_jd(1: 4);

                                % negate kd flag at those 4 pts AND store those 4
                                % pts in correspondence to i j k point for later
                                % use
                                k_interp_point_count = k_interp_point_count + 1;
                                k_interp_points{k_interp_point_count} = [i j];
                                source_ids = zeros(4, 2);
                                for l = 1: 4
                                    obj.grids{kd}.flag(closest_id(l), closest_jd(l)) = -abs(obj.grids{kd}.flag(closest_id(l), closest_jd(l)));
                                    source_ids(l, :) = [closest_id(l) closest_jd(l)];
                                end
                                k_interp_source_ids{k_interp_point_count} = source_ids;

                                obj.grids{k}.flag(i, j) = -abs(obj.grids{k}.flag(i, j));
                                % mark points on upper grids needed for interpolation by lower grid
                                %for id = 1: 4
                                %    for jd = 1: 4
                                %        %grids{kd}.flag(closest_id(id), closest_jd(jd)) = -k;
                                %    end
                                %end
                            else
                                obj.grids{k}.flag(i, j) = 0;
                            end
                        end
                    end
                end

                obj.grids{k}.interp_point_count = obj.grids{k}.interp_point_count + k_interp_point_count;
                obj.grids{k}.interp_points = [obj.grids{k}.interp_points, k_interp_points];
                obj.grids{k}.interp_source_ids = [obj.grids{k}.interp_source_ids, k_interp_source_ids];

                % add code for "if grids{k}.flag(i, j) can be interior point"

            end

            % STEP 6: Correct signs
            for k = 1: obj.n_grids
                for i = 1: obj.grids{k}.ny
                    for j = 1: obj.grids{k}.nx
                        if abs(obj.grids{k}.flag(i, j)) == k % discretization point
                            obj.grids{k}.flag(i, j) = abs(obj.grids{k}.flag(i, j));
                        elseif abs(obj.grids{k}.flag(i, j)) > 0 % interpolation point
                            obj.grids{k}.flag(i, j) = -abs(obj.grids{k}.flag(i, j));
                        end
                    end
                end
            end
        end
        
        % compute interpolated values at grids
        function [] = interpolate(obj)
            for k = 1: obj.n_grids
                disp(strcat('overset: interpolating at grid ', obj.grids{k}.name));
                for i = 1: obj.grids{k}.interp_point_count
                    interp_point_ids = obj.grids{k}.interp_points(i);
                    interp_kd = -1 * obj.grids{k}.flag(interp_point_ids{:}(1), interp_point_ids{:}(2));
                    interp_point = obj.grids{k}.get_global_coords_at(interp_point_ids{:}(1), interp_point_ids{:}(2));
                    interp_source_ids = obj.grids{k}.interp_source_ids(i);
                    
                    source_coords = zeros(4, 2);
                    source_vals = zeros(4, 1);

                    for j = 1: 4
                        source_coords(j, :) = obj.grids{interp_kd}.get_global_coords_at(interp_source_ids{:}(j, 1), interp_source_ids{:}(j, 2));
                        source_vals(j) = obj.grids{interp_kd}.val(interp_source_ids{:}(j, 1), interp_source_ids{:}(j, 2));
                    end
                    
                    interp_val = griddata(source_coords(:, 2), source_coords(:, 1), source_vals, interp_point(2), interp_point(1), 'linear');
                    obj.grids{k}.val(interp_point_ids{:}(1), interp_point_ids{:}(2)) = interp_val;
                end
            end
        end
        
        % display grid
        function [] = display_grid(obj, fig)
            disp(strcat('overset: printing composite grid ', obj.name));
            for k = 1: obj.n_grids
                obj.grids{k}.display_grid(fig);
            end
        end
        
        % display data
        function [] = display_data(obj, fig)
            disp(strcat('overset: printing data on composite grid ', obj.name));
            for k = 1: obj.n_grids
                obj.grids{k}.display_data(fig);
            end
        end
        
    end
end