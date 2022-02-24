%% This code grids the depth of IRHs over the WAIS and exports as tif file
%
% Written by J. Bodart (UoE) - last updated: 15/02/2021
%
%%
clear all
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));

%% load and extract IRHs boundary shapefile
% shapefile = shaperead('D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\Outlines_IRHs\PIG_IMAFI_IRH_boundary.shp');
% shp_x = shapefile.X;
% shp_y = shapefile.Y;

%% import R2 Bodart
% PASIN
H2_Bodart=['D:\University_Edinburgh\PhD\PIG_manuscript\PIG_review\picked_layers\Individual_updated_EPSG3031\Bodart_IRHs_dataset_final\R2_PASIN.csv'];
fid=fopen(H2_Bodart);
H2_Bodart=textscan(fid,'%f %f %f %f %f %f %f %f','delimiter',',','headerLines', 2);
fclose(fid);

R2_x =(H2_Bodart{1,1});
R2_y =(H2_Bodart{1,2});
R2_depth =(H2_Bodart{1,6});

% OIB 1
H2_Bodart=['D:\University_Edinburgh\PhD\PIG_manuscript\PIG_review\picked_layers\Individual_updated_EPSG3031\Bodart_IRHs_dataset_final\R2_OIB_20181022_01.csv'];
fid=fopen(H2_Bodart);
H2_Bodart=textscan(fid,'%f %f %f %f %f %f %f %s','delimiter',',','headerLines', 2);
fclose(fid);

R2_x = vertcat(R2_x, H2_Bodart{1,1});
R2_y = vertcat(R2_y, H2_Bodart{1,2});
R2_depth = vertcat(R2_depth, H2_Bodart{1,6});

% OIB 2
H2_Bodart=['D:\University_Edinburgh\PhD\PIG_manuscript\PIG_review\picked_layers\Individual_updated_EPSG3031\Bodart_IRHs_dataset_final\R2_OIB_20161109_03.csv'];
fid=fopen(H2_Bodart);
H2_Bodart=textscan(fid,'%f %f %f %f %f %f %f %s','delimiter',',','headerLines', 2);
fclose(fid);

R2_x = vertcat(R2_x, H2_Bodart{1,1});
R2_y = vertcat(R2_y, H2_Bodart{1,2});
R2_depth = vertcat(R2_depth, H2_Bodart{1,6});

%% import H2 Ashmore
H2_Ashmore=['D:\University_Edinburgh\AntArchitecture\Ashmore_2020_data\v2\H2_Raw_txt.txt'];
fid=fopen(H2_Ashmore);
H2_Ashmore=textscan(fid,'%f %f %f %f %f %s','delimiter',' ','headerLines', 1);
fclose(fid);

H2_x =(H2_Ashmore{1,2});
H2_y =(H2_Ashmore{1,3});
H2_depth =(H2_Ashmore{1,5});

%% combine Bodart and Ashmore arrays
R2_x = vertcat(R2_x,H2_x);
R2_y = vertcat(R2_y,H2_y);
R2_depth = vertcat(R2_depth,H2_depth);

%% import LM9 Muldoon
H2_Muldoon = ['D:\University_Texas\UTIG_layers\UTIG_layers\lm-MERGE-lay9-grg-depth_EPSG3031.txt'];
fid=fopen(H2_Muldoon);
H2_Muldoon=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines', 1);
fclose(fid);

LM9_x =(H2_Muldoon{1,2});
LM9_y =(H2_Muldoon{1,3});
LM9_depth =(H2_Muldoon{1,4});

%% import Shallow 2 Leuro
% S2_Leuro=['D:\University_Texas\UTIG_layers\UTIG_layers\Leuro_layers\THW_LEURO_SHALLOW_2.csv'];
% fid=fopen(S2_Leuro);
% S2_Leuro=textscan(fid,'%s %f %f %f %f %f %f %f %f %s %s','delimiter',',','headerLines', 24);
% fclose(fid);
% 
% R2_x =(S2_Leuro{1,4});
% R2_y =(S2_Leuro{1,5});
% R2_depth =(S2_Leuro{1,9});
% R2_depth (R2_depth < 0) = NaN; % for Thwaites Leuro to exclude outliers

%% Combine Bodart, Ashmore and Muldoon IRHs
R2_x = vertcat(R2_x,LM9_x);
R2_y = vertcat(R2_y,LM9_y);
R2_depth = vertcat(R2_depth,LM9_depth);

%% combine arrays
R2_array = [R2_x, R2_y, R2_depth];

% remove nan
nans = all(all(isnan(R2_array),3),2);
R2_array(nans,:,:) = [];

%% create grid with 1-km spacing around West Antarctica
%[x,y] = psgrid(-77.34,-95.56,600,1,'xy'); % Centered over PIG - resolution grid: 1km
%[x,y] = psgrid(-77.5820,-108.4458,[250 420],5,'xy'); % Centered vertically to get PIG and THW Bodart and Muldoon IRHs - resolution is 5 km grid
%[x,y] = psgrid(-79.5820,-93.56,[1300 750],1,'xy'); % Centered horizontally to get PIG and IMAFI Bodart and Ashmore IRHs - resolution is 1 km grid
%[x,y] = psgrid(-77.3807,-112.4346,[650 750],1,'xy'); % Centered around Thwaites for Leuro layers - resolution is 1 km grid
[x,y] = psgrid(-1.16452e+06,-3.07807e+05,[1100 1200],1,'xy'); % Centered vertically to get Bodart, Ashmore and Muldoon IRHs - resolution is 1 km grid

% figure
% bedmap2('patchshelves','xy')
% hold on
% bedmap2('patchgl','xy')
% plot(x,y,'b.')
% hold on
% plot(shp_x,shp_y)
% plot(R2_x,R2_y)

%% create interpolation from scatter with original XY data
%R2_interp = scatteredInterpolant(R2_x,R2_y,R2_depth);

% fit interpolated data to grid
%R2_grid = R2_interp(x,y);
R2_grid = griddata(R2_array(:,1),R2_array(:,2),R2_array(:,3),x,y);

% smooth data using gaussian filter
R2_grid_smooth = smoothdata(R2_grid,'gaussian',20);

%% plot
% set axes
%axes = ([-1.801e6 -0.91e6 -5e5 2e5]); % PIG
axes = ([-1.6e6 -0.6e6 -9.5e5 2.5e5]); % PIG - IMAFI
%axes = ([-1.65e6 -0.9e6 -9.5e5 -1.5e5]); % Thwaites
% 
% % plot figure
figure;
axis (axes) % set figure axes
bedmap2('patchshelves','xy')
hold on
plot(x,y,'b.')
imagescn(x,y,R2_grid_smooth);
hold on
plot(shp_x,shp_y)
scatter(R2_x,R2_y,1,'k','o','filled')
scatter(H2_x,H2_y,1,'r','o','filled')
colormap(flipud(parula))
cb = colorbar;
set(cb, 'YDir', 'reverse' );
caxis([200 2200])
title('Natural - gaussian 20')
% plot bedmap bed for PIG 
% figure;
% axis (axes) % set figure axes
% bedmap2('bed','parula','xy') 
% hold on
% scatter(R2_x,R2_y,1,'k','o','filled')
% scatter(H2_x,H2_y,1,'r','o','filled')

%% export grid as geotiff file
cd 'D:\R_University_Edinburgh\WAIS_accumulation\grid_IRHs'
%[H2_Bodart,R] = geotiffread('D:\R_University_Edinburgh\WAIS_accumulation\H2_raster_points\H2_point_raster.tif');

% get xy limits from gridded map
xlimits = [min(x(1,:)) max(x(1,:))];
ylimits = [min(y(1,:)) max(y(end,:))];

% specify R object
R = maprefcells(xlimits,ylimits,size(R2_grid_smooth), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');

coordRefSysCode = 3031; % EPSG
filename = 'R2_gridded_PIG_IMAFI_THW_smooth_1km.tif';
geotiffwrite(filename,R2_grid_smooth,R,'CoordRefSysCode', coordRefSysCode);

%%