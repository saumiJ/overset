classdef overset_base_grid < overset_grid
% class for the overset base-grid

    properties
        % see overset_grid.m
    end
    
    methods
        
        function obj = overset_base_grid(name_, id_, ...
                nx_, ny_, ...
                dx_, dy_, ...
                sg1_list_voidBoundaryPointList_) % constructor
            % call base-class constructor
            obj@overset_grid(name_, id_, nx_, ny_, dx_, dy_, sg1_list_voidBoundaryPointList_); 
        end
        
        function global_coords = get_global_coords(obj) % returns global coordinates of grid points
            global_coords = obj.grid_coords;
        end
        
        function [poly_x, poly_y] = get_void_polygon(obj, k)
            poly_x = zeros(1, size(obj.void_polygons{k}, 2));
            poly_y = zeros(1, size(obj.void_polygons{k}, 2));
            for l = 1: size(obj.void_polygons{k}, 2)
                poly_x(l) = obj.dx*(obj.void_polygons{k}(2, l) - 1);
                poly_y(l) = obj.dy*(obj.void_polygons{k}(1, l) - 1);
            end
        end
        
        function fig = display_grid(obj, fig) % plots the grid on figure img_number
            figure(fig);
            hold on
            y = zeros(1, obj.nx*obj.ny);
            x = zeros(1, obj.nx*obj.ny);

            k = 1;

            for i = 1: obj.ny
                for j = 1: obj.nx
                    if ~obj.isVoidBoundary(i, j) && obj.flag(i, j) ~= 0
                        y(1, k) = obj.grid_coords(i, j, 1);
                        x(1, k) = obj.grid_coords(i, j, 2);
                        k = k + 1;
                    end
                end
            end
            
            scatter(x, y, 12, 'filled');
            hold off;
        end
        
    end
end