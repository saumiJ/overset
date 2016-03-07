classdef overset_grid < handle
% parent class for overset grids
    properties

        % TODO:
        % - convert representation of points to a class
        
        % general properties
        name;   % data name
        id;     % grid id
        nx;     % grid points in x
        ny;     % grid points in y
        dx;     % grid spacing in x
        dy;     % grid spacing in y

        % data properties
        val;        % data
        flag;  % is grid point active?
        isVoidBoundary;
        num_void_polygons;
        void_polygons;

        % spatial properties
        grid_center;        % mid-point of grid
        grid_coords;        % local coords of grid
        
        % interpolation points details
        interp_point_count
        interp_points % indices
        interp_source_ids
    
    end
    
    methods
        % constructor
        function obj = overset_grid(name_, ...
                id_, ...
                nx_, ny_, ...
                dx_, dy_ ,...
                list_voidBoundaryPointList_ ...
                )
            obj.name = name_;
            obj.id = id_;
            obj.nx = nx_; obj.ny = ny_;
            obj.dx = dx_; obj.dy = dy_;
            
            obj.val = zeros(obj.ny, obj.nx); 
            obj.flag = -0.001*ones(obj.ny, obj.nx); % assign an invalid value
            obj.isVoidBoundary = zeros(obj.ny, obj.nx); % no voids by default
            obj.grid_center = [((obj.ny-1) * obj.dx)/2 ((obj.nx-1) * obj.dy)/2];
            obj.grid_coords = zeros(obj.ny, obj.nx, 2);
            
            for i = 1: obj.ny
                for j = 1: obj.nx
                    obj.val(i, j) = i+j;
                    obj.grid_coords(i, j, 1) = (i-1)*obj.dy;
                    obj.grid_coords(i, j, 2) = (j-1)*obj.dx;
                end
            end
            
            obj.num_void_polygons = length(list_voidBoundaryPointList_);
            obj.void_polygons = cell(length(list_voidBoundaryPointList_));
            for k = 1: obj.num_void_polygons
                obj.void_polygons{k} = list_voidBoundaryPointList_{k};
                for l = 1: size(list_voidBoundaryPointList_{k}, 2)
                    obj.isVoidBoundary(list_voidBoundaryPointList_{k}(1, l), ...
                        list_voidBoundaryPointList_{k}(2, l)) = true;
                end
            end
            
        end
    end
end