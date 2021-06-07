function [grid_voronoi]=voronoi_func_v4(y_coor,x_coor,data,gridY,gridX,Epsilon)
% INPUTS:
%   y_coor: network stations latitude (WGS)  [vector]
%   x_coor: network stations longitude (WGS) [vector]
%   data: individual station SWC data [vector of (x_coor,y_coor)]
%   gridY: latitude of large area needed for initial voronoi
%   gridX: longitude of larger area needed for initial voronoi
%   Epsilon: 1(UTM) or 1E-8 (WGS)

% OUTPUTS
%   grid_voronoi:raster of the station index assigned to each grid cell 


%Find the stations in the larger grid. The SquareBV function does not produce good results if there is points outside the grid
ind = find(x_coor>min(min(gridX)) & x_coor<max(max(gridX)) & y_coor>min(min(gridY)) & y_coor<max(max(gridY)) & ~isnan(data));

station_number=[1:length(data)]';

if length(ind) < 3 % if there are only 2 stations then this function does work
    data_out     = NaN;
    data_out_std = NaN;    
else
    %reduce station longitude, latitude, and data to those within the larger area needed for the initial voronoi
    x_coor = x_coor(ind); y_coor=y_coor(ind); data=data(ind); station_number=station_number(ind);
    
    %initial voronoi area
    [~,vpoly]=SquareBV(x_coor,y_coor,0,Epsilon,[min(min(gridX)),max(max(gridX)),min(min(gridY)),max(max(gridY))]);
    %cell area and polygons of the cell(vpoly) are the outputs of SquareBV
      
    % Initialize raster of voronoi and convert ouputed polygons to raster
    grid_voronoi = gridX*0;
    for i=1:length(vpoly)
        %find the raster grid points within each voronoi polygon and label them
        in = inpolygon(gridX,gridY,vpoly{1,i}(:,1),vpoly{1,i}(:,2)); hold on
        grid_voronoi(in)=station_number(i);
    end
    
end
