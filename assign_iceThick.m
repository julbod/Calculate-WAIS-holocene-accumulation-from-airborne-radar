%% This code is an alternative to the more elaborate 'NNjoin_assign_iceThick'
% MATLAB code. Instead of doing the NNjoin operation in the code, this code
% uses data that has already been joined in QGIS using NNjoin tool. This
% data is imported here and ice thickness values from BedMachine is
% assigned to ice thickness values which are either 450-m away from the
% nearest ice thickness point from the radar flightline (e.g. IMAFI), 
% or that is NaN as a result of the missing thickness value from the radar 
% data (e.g. BBAS). 
%
% The code reads in the tabular data for the IRH, if a 'distance' variable
% exists in the file/code (e.g. IMAFI), the code places '-9999' for values
% which are 450-m away from nearest ice thickness point, or skips this
% step if this has already been done previously and contains nan where
% there is no radar ice thickness (e.g. BBAS). The code then imports the
% BedMachine tif file, extracts the xyz data, and then for all -9999 or NaN
% values finds the nearest BedMachine node and assign the ice thickness
% value to these points. The code then exports the data as text file. 
%
% Written by J. Bodart (UoE) - 23.02.2022
%
%%
clear all
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));
addpath 'D:\R_University_Edinburgh\NKK_DWA_JAB_correlation'

%% import IRH
layer=['D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Ashmore_IRHs\H3_v2_iceThick.csv']; % Ashmore
%layer=['D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Muldoon_IRHs\LM9_EPSG3031.csv']; % Muldoon
%layer=['D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Bodart_IRHs\R3_Bodart_all.txt']; % Bodart
fid=fopen(layer);
layer=textscan(fid,'%f %f %f %f %s %f %f','delimiter',',','headerLines',1); % Ashmore
%layer=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1); % Muldoon
%layer=textscan(fid,'%f %f %f %f','delimiter',',','headerLines',1); % Bodart
fclose(fid);

layer_x =(layer{1,1});
layer_y =(layer{1,2});
layer_depth =(layer{1,3});
layer_iceThick =(layer{1,4});
layer_dist =(layer{1,7});

% replace nan by -9999 if any
layer_iceThick(isnan(layer_iceThick)) = -9999;

% check distribution of data using own function
quart_layer = Quartiles_funcs(layer_dist);

% concat xy data
layer_xy = horzcat(layer_x,layer_y);

%% Place -9999 if distance is larger than 450 m
if exist('layer_dist','var')
    for i = 1:length(layer_dist)
        i
        % make sure minimum distance between points is less than 450
        if layer_dist(i) > 500
            layer_iceThick (i) = -9999; % if larger than 450 m, place -9999
        end
    end
end

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
thick = geotiffread(tif);

% flip array in up/down direction
thick = flipud(thick);

% conver to vector
grid_data = [x_grd(:) y_grd(:) thick(:)];
clear thick

%% Assign BedMachine data points to all -9999 values
for i = 1:length(layer_iceThick)
    if layer_iceThick (i) == -9999
        % compute eucledian distance
        dist = sqrt(sum(bsxfun(@minus, grid_data(:,1:2), layer_xy(i,1:2)).^2,2));

        % find closest value/index from distance
        %clos_idx (i) = find(dist==min(dist));
        
        % get matching ice thickness from  BedMachine
        thick (i) = grid_data(find(dist==min(dist),1),3);
    end
    
    if i == round(length(layer_iceThick)*0.25)
        display('#### 1/4 of the way there #### ')
    elseif i == round(length(layer_iceThick)*0.5)
        display('#### 2/4 of the way there #### ')
    elseif i == round(length(layer_iceThick)*0.75)
        display('#### 3/4 of the way there #### ')
    end
end

%% place data into original array
layer_iceThick(layer_iceThick == -9999) = thick(find(thick ~= 0));

%% export data to csv
if exist('layer_dist','var')
    table = table(layer_xy(:,1),layer_xy(:,2),layer_depth, layer_iceThick, layer_dist, 'VariableNames', { 'x', 'y','lyr_depth','iceThick', 'distance'} );
    writetable(table, 'D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Ashmore_IRHs\H3_iceThick_final.txt')
else
    table = table(layer_xy(:,1),layer_xy(:,2),layer_depth, layer_iceThick, 'VariableNames', { 'x', 'y','lyr_depth','iceThick'} );
    writetable(table, 'D:\R_University_Edinburgh\WAIS_accumulation\IRHs\Ashmore_IRHs\H3_iceThick_final.txt')
