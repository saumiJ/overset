classdef overset_base_grid < overset_grid
% class for the overset base-grid

    properties
        % see overset_grid.m
    end
    
    methods
        
        function obj = overset_base_grid(name_, ...
                nx_, ny_, ...
                dx_, dy_) % constructor
            % call base-class constructor
            obj@overset_grid(name_, nx_, ny_, dx_, dy_); 
        end
        
        function fig = display_grid(obj, fig) % plots the grid on figure img_number
            figure(fig);
            hold on
            y = zeros(1, obj.nx*obj.ny);
            x = zeros(1, obj.nx*obj.ny);

            k = 1;

            for i = 1: obj.ny
                for j = 1: obj.nx
                    y(1, k) = obj.grid_coords(i, j, 1);
                    x(1, k) = obj.grid_coords(i, j, 2);
                    k = k + 1;
                end
            end
            
            scatter(x, y, 12, 'filled');
            hold off;
        end
        
    end
end