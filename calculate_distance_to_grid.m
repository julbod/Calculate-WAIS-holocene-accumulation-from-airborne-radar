%% This code calculates the distance between a grid cell in the gridded
% IRH file from Bodart et al. and the nearest IRH data point along a 
% radar track line.
%
% Code last updated 01/03/2022 by J. Bodart (UoE)
%
%% Code explanation
% The code first imports the gridded IRH accumulation tif file from Bodart 
% et al., extracts the xy information, and then imports the IRH files for
% R2 over IMAFI, PIG and THW. The code then find the minimum distance
% between the grid cell and the nearest IRH along a radar track line. This
% allows users to know how far the interpolation is going and where it is
% likely to be least accurate.
%
%% set variables and directories
clear all
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));

%% load grid
% load gridded IRH accumulation data from this study
% Reference: Bodart et al., 2022
% Details: 4.72 ka; units: m/yr ice equivalent; grid res: 1 km
% Download: 
tif = 'D:\R_University_Edinburgh\WAIS_accumulation\aligned_grids\R2_PIG-IMAFI-THW_ACCU_NYE_1_GRDSMOOTH_ALIGNED.tif';
tif_info = geotiffinfo(tif);

% extract x and y coordinates from tif
[x_min, y_min] = deal(tif_info.BoundingBox(1, 1), tif_info.BoundingBox(1, 2));
[num_x, num_y] = deal(tif_info.Width, tif_info.Height);
tif_inc = tif_info.GeoTIFFTags.ModelPixelScaleTag(1);

% assign coordinate data
x_coords = (x_min + (tif_inc / 2)) + (0:tif_inc:((num_x - 1) * tif_inc));
y_coords = (y_min + (tif_inc / 2)) + (0:tif_inc:((num_y - 1) * tif_inc))';

% grid xy coordinates in meters
[x_grd, y_grd] = meshgrid(x_coords, y_coords);

% make column vectors
x_grd = x_grd(:);
y_grd = y_grd(:);
xy_grd = horzcat(x_grd, y_grd);

% get data
[grid_vals, R_grid, grids] = geotiffread(tif);
grid_vals (grid_vals==0) = NaN;
grid_vals=double(grid_vals);

%% load IRHs 
IMAFI = ['D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Ashmore_IRHs\H2_iceThick_final.txt'];
fid=fopen(IMAFI);
IMAFI=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);

IMAFI_x =(IMAFI{1,1});
IMAFI_y =(IMAFI{1,2});
IMAFI_depth =(IMAFI{1,3});
IMAFI_iceThick =(IMAFI{1,4});
IMAFI = horzcat(IMAFI_x,IMAFI_y,IMAFI_depth,IMAFI_iceThick);

PIG =['D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Bodart_IRHs\R2_all_iceThick_final.txt'];
fid=fopen(PIG);
PIG=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);

PIG_x =(PIG{1,1});
PIG_y =(PIG{1,2});
PIG_depth =(PIG{1,3});
PIG_iceThick =(PIG{1,4});
PIG = horzcat(PIG_x,PIG_y,PIG_depth,PIG_iceThick);

THW =['D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Muldoon_IRHs\LM9_iceThick_final.txt'];
fid=fopen(THW);
THW=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);

THW_x =(THW{1,1});
THW_y =(THW{1,2});
THW_depth =(THW{1,3});
THW_iceThick =(THW{1,4});
THW = horzcat(THW_x,THW_y,THW_depth,THW_iceThick);

%% join layers together
R2_array = vertcat(IMAFI,PIG,THW);

%% calculate distance
disp('Calculating minimum distance for each grid cell...')

%f = waitbar(0,'Starting...');
%for i = 1:length(xy_grd)
for i = 1:length(xy_grd)
    
    % get 100% waitbar
    %waitbar(i/length(layer_xy), f, sprintf('Progress: %d %%', floor(i/length(layer_xy)*100)));
    %pause(0.1);
    %i
    % compute eucledian distance
    distances = sqrt(sum(bsxfun(@minus, R2_array(:,1:2), xy_grd(i,1:2)).^2,2));
    
    % find closest value/index from distance
    %closest_value = thick_xy(find(distances==min(distances)),:);
    closest_idx (i) = find(distances==min(distances),1); 
    distance (i) = min(distances);
    
    if i == round(length(xy_grd)*0.25)
        display('#### 1/4 of the way there #### ')
    elseif i == round(length(xy_grd)*0.5)
        display('#### 2/4 of the way there #### ')
    elseif i == round(length(xy_grd)*0.75)
        display('#### 3/4 of the way there #### ')
    end
end

distance = distance./1e4; % convert from meters to km

%% create grid with 1-km spacing around West Antarctica
%[x,y] = psgrid(-801500, -301500,[2084 1406],1,'xy'); % Centered vertically to get Bodart, Ashmore and Muldoon IRHs - resolution is 1 km grid
% alternative to above:
% x = grids(1,1) : 1000 : grids(1,2); 
% y = grids(2,1) : 1000 : grids(2,2); 
% [x,y] = meshgrid(x,y); 

%% create interpolation from scatter with original XY data
% fit interpolated data to grid
%dist_grid = griddata(xy_grd(:,1),xy_grd(:,2),distss,x,y);

%% write data to GeoTIF
% get xy limits from gridded map
%cd 'D:\R_University_Edinburgh\WAIS_accumulation\grid_IRHs_depth\calculate_distance'
% xlimits = [min(x(1,:)) max(x(1,:))];
% ylimits = [min(y(1,:)) max(y(end,:))];
% 
% % specify R object
% R = maprefcells(xlimits,ylimits,size(dist_grid), ...
%     'ColumnsStartFrom','south','RowsStartFrom','west');
% 
% coordRefSysCode = 3031; % EPSG
% filename = 'R2_grid_distance_km.tif'; 
% geotiffwrite(filename, dist_grid,R,'CoordRefSysCode', coordRefSysCode);

%% export line data to csv
table = table(xy_grd(:,1),xy_grd(:,2),distance, 'VariableNames', { 'x', 'y','distance_km'} );
writetable(table, 'D:\R_University_Edinburgh\WAIS_accumulation\grid_IRHs_depth\calculate_distance\R2_grid_distance_km.txt')

%%