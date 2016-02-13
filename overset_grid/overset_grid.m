classdef overset_grid
% parent class for overset grids
    properties

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

        % spatial properties
        grid_center;        % mid-point of grid
        grid_coords;        % local coords of grid
    
    end
    
    methods
        function obj = overset_grid(name_, ...
                id_, ...
                nx_, ny_, ...
                dx_, dy_ ...
                ) % constructor
            obj.name = name_;
            obj.id = id_;
            obj.nx = nx_; obj.ny = ny_;
            obj.dx = dx_; obj.dy = dy_;
            
            obj.val = zeros(obj.ny, obj.nx); 
            obj.flag = -ones(obj.ny, obj.nx);
            obj.grid_center = [((obj.ny-1) * obj.dx)/2 ((obj.nx-1) * obj.dy)/2];
            obj.grid_coords = zeros(obj.ny, obj.nx, 2);
            for i = 1: obj.ny
                for j = 1: obj.nx
                    obj.grid_coords(i, j, 1) = (i-1)*obj.dy;
                    obj.grid_coords(i, j, 2) = (j-1)*obj.dx;
                end
            end
            
        end
    end
end