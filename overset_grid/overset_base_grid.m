classdef overset_base_grid < overset_grid
% class for the overset base-grid

    properties
        % see overset_grid.m
    end
    
    methods
        
        % constructor
        function obj = overset_base_grid(name_, id_, ...
                nx_, ny_, ...
                dx_, dy_, ...
                sg1_list_voidBoundaryPointList_)
            disp(strcat('overset: constructing base grid ', name_));
            % call base-class constructor
            obj@overset_grid(name_, id_, nx_, ny_, dx_, dy_, sg1_list_voidBoundaryPointList_); 
        end
        
        % get global coordinates of a point
        function global_coords_at = get_global_coords_at(obj, i, j)
            global_coords_at = zeros(1, 2);
            
            global_coords_at(1, 1) = obj.grid_coords(i, j, 1);
            global_coords_at(1, 2) = obj.grid_coords(i, j, 2);
        end
        
        % get global coordinates
        function global_coords = get_global_coords(obj) % returns global coordinates of grid points
            global_coords = obj.grid_coords;
        end
        
        % get k-th void polygon
        function [poly_x, poly_y] = get_void_polygon(obj, k)
            poly_x = zeros(1, size(obj.void_polygons{k}, 2));
            poly_y = zeros(1, size(obj.void_polygons{k}, 2));
            for l = 1: size(obj.void_polygons{k}, 2)
                poly_x(l) = obj.dx*(obj.void_polygons{k}(2, l) - 1);
                poly_y(l) = obj.dy*(obj.void_polygons{k}(1, l) - 1);
            end
        end
        
        % plots the grid on figure fig
        function fig = display_grid(obj, fig) 
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
            
            disp(strcat('overset: printing base grid ', obj.name));
            scatter(x, y, 12, 'filled');
            hold off;
        end
        
        % displays data as a contour on figure fig
        function [] = display_data(obj, fig)
            figure(fig);
            hold on
            y = zeros(1, obj.ny * obj.nx);
            x = zeros(1, obj.ny * obj.nx);
            data = zeros(1, obj.ny * obj.nx);
            
            k = 1;
            for i = 1: obj.ny
                for j = 1: obj.nx
                    if ~obj.isVoidBoundary(i, j) && obj.flag(i, j) == obj.id+1
                    %if ~obj.isVoidBoundary(i, j) && obj.flag(i, j) ~= 0
                        y(1, k) = obj.grid_coords(i, j, 1);
                        x(1, k) = obj.grid_coords(i, j, 2);
                        data(1, k) = obj.val(i, j);
                        k = k + 1;
                    end
                end
            end
            
            disp(strcat('overset: printing data on base grid ', obj.name));
            scatter3(x, y, data, 20, data, 'filled');
            
            hold off
        end
    end
end