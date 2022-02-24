%% This code joins multiple datasets together based on spatial location
% In effect, this code is an alternative to the QGIS plugin NNJoin
% <https://plugins.qgis.org/plugins/NNJoin/> but the key difference to this
% plugin is that here we are able to set a maximum distance between points
% beyond which the spatial join is not made an NaNs are assigned. For those
% values which are outside of the desired distance, the code then (step 2)
% imports the BedMachine ice thickness products and assigns the ice thickness
% value for all points falling within the nearest gridded-product node. A
% simple alternative to this code, which does step 2 but skips step 1 can
% be found as 'assign_iceThick.m'. 
%
% First, the code imports the IRH depth data and the radar ice thickness
% point data and calculates the distance between all the points in the radar 
% ice thickness data and each point in the IRH depth data. The minimum
% distance is then found, and the resulting ice thickness from the nearest
% point is calculated. If the minimum distance is greater than the cell
% resolution of the BedMachine gridded product (~450 m), the code assigns
% NaNs (or -9999). Following this, the BedMachine ice thickness product
% is imported and ice thickness values from BedMachine are assigned to ice
% thickness values which are either 450-m away from the nearest ice 
% thickness point from the radar flightline. Finally, the code exports the 
% data as tabular data file. 
%
% Written by J. Bodart (UoE) - 23.02.2022
%
%%
clear all
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));

%% import IRH
layer=['D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Muldoon_IRHs\LM9_EPSG3031.csv']; % Muldoon
fid=fopen(layer);
layer=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1); %Muldoon
fclose(fid);

layer_x =(layer{1,2});
layer_y =(layer{1,3});
layer_depth =(layer{1,4});

%% import ice thickness vector
thick=['D:\University_Texas\UTIG_YYYY_AGASEA-BAS_AIR_BM2.csv'];
fid=fopen(thick);
thick=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %f','delimiter',',','headerLines', 1); % Muldoon
fclose(fid);

thick_lon =(thick{1,1});
thick_lat =(thick{1,2});
thick_thick =(thick{1,5});

% convert lon/lat to psx/y
[thick_x thick_y] = ll2ps(thick_lat, thick_lon);

clear thick layer

%% assign ice thickness values to closest point in IRHs
layer_xy = horzcat(layer_x,layer_y);
thick_xy = horzcat(thick_x,thick_y);

%f = waitbar(0,'Starting...');
for i = 1:length(layer_xy)
    
    % get 100% waitbar
    %waitbar(i/length(layer_xy), f, sprintf('Progress: %d %%', floor(i/length(layer_xy)*100)));
    %pause(0.1);
    %i
    % compute eucledian distance
    distances = sqrt(sum(bsxfun(@minus, thick_xy(:,1:2), layer_xy(i,1:2)).^2,2));
    
    % find closest value/index from distance
    %closest_value = thick_xy(find(distances==min(distances)),:);
    closest_idx (i) = find(distances==min(distances));
    distance (i) = min(distances);
    
    % make sure minimum distance between points is less than 450
    if min(distances) > 450
        thick (i) = -9999; % if larger than 450 m, place -9999
    else
        thick (i) = thick_thick(closest_idx(i)); % if smaller, get thickness value
    end
    
    if i == round(length(layer_xy)*0.25)
        display('#### 1/4 of the way there #### ')
    elseif i == round(length(layer_xy)*0.5)
        display('#### 2/4 of the way there #### ')
    elseif i == round(length(layer_xy)*0.75)
        display('#### 3/4 of the way there #### ')
    end
    
end

thick_radar = thick;

%% read BedMachine tif file data
tif = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\BedMachine_thick_clipped_450m.tif';
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
[num_y, num_x] = size(x_grd);

% read thickness data
grd_thick = geotiffread(tif);

% flip array in up/down direction
BedMac_thick = flipud(grd_thick);

% conver to vector
grid_data = [x_grd(:) y_grd(:) grd_thick(:)];

%%
%% Assign BedMachine data points to all -9999 values
for i = 1 :length(thick)
    if thick (i) == -9999
        
        % compute eucledian distance
        dist = sqrt(sum(bsxfun(@minus, grid_data(:,1:2), layer_xy(i,1:2)).^2,2));

        % find closest value/index from distance
        %clos_idx (i) = find(dist==min(dist));
        
        % get matching ice thickness from  BedMachine
        thick (i) = grid_data(find(dist==min(dist),1),3);
    end
    
    if i == round(length(thick)*0.25)
        display('#### 1/4 of the way there #### ')
    elseif i == round(length(thick)*0.5)
        display('#### 2/4 of the way there #### ')
    elseif i == round(length(thick)*0.75)
        display('#### 3/4 of the way there #### ')
    end
end

%% export data to csv
table = table(layer_xy(:,1),layer_xy(:,2),layer_depth, thick, distance.', 'VariableNames', { 'x', 'y','lyr_depth','iceThick','distance'} );
writetable(table, 'D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Muldoon_IRHs\LM9_iceThick_final.txt')

%%