%% This code calculates vertical strain rates from Nye-calculated accumulation
% for the 4.72 ka IRH over the WAIS see Bodart et al. (2023, The Cryosphere)
% for more details.
%
% Explanation of code:
% First, the code imports the text file containing the accumulation rates
% calculated using 'Calculate_accumulation_rates.m' file for the 4.72 ka IRH
% over the WAIS, and calculates vertical strain rates using a simple 
% equation. The code finally grids and smoothes the data using gaussian 
% filtering and finally grids and exports the data as 
% geotiff files for import into QGIS and to plot as figures. This data is
% shown in Fig. S2a of Bodart et al. (2023).
% 
% Note that some functions used here originate from the Antarctic Mapping
% Tools toolbox developped by Chad Green (and available at:
% https://github.com/chadagreene/Antarctic-Mapping-Tools)
%
% Code written by J. Bodart (UoE) - 23/02/2022
%
%%
clear all

%% import layer data
R2_accu = ['...\accumulation_IRH_combined.txt'];
fid=fopen(R2_accu);
R2_accu=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);

x_accu =(R2_accu{1,1}); % PSX
y_accu =(R2_accu{1,2}); % PSY
z_depth =(R2_accu{1,3}); % Depth of IRH
iceThick =(R2_accu{1,4}); % Ice thickness at IRH point
accu =(R2_accu{1,5}); % Nye accumulation

% Combine into one array
R2_accu = horzcat(x_accu,y_accu,z_depth,iceThick,accu);

%% Calculate vertical strain rates using Eq. 2 in paper
strain = (-(accu))./iceThick;

%% Create grid with 1-km spacing around West Antarctica
% centered vertically to get all WAIS IRH coverage
[x,y] = psgrid(-801500, -301500,[2084 1406],1,'xy');

%% Apply filter and grid the along-track data
% smooth along-track variables using moving-average gaussian filter
strain_smooth = smoothdata(strain,'gaussian',30);

% interpolate filtered along-track data to grid
strain_grid = griddata(R2_accu(:,1),R2_accu(:,2),strain_smooth,x,y, 'natural');

% filter grid using average cell filter
filter = fspecial('average',[18,18]);
strain_grid_smooth = imfilter(strain_grid,filter);

%% export filtered grid data as geotiff file for QGIS import and figures
% get xy limits from grid
xlimits = [min(x(1,:)) max(x(1,:))];
ylimits = [min(y(1,:)) max(y(end,:))];

% specify R object for Geotiff
R = maprefcells(xlimits,ylimits,size(strain_grid_smooth), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');

% specify coordinate system and name of files
coordRefSysCode = 3031; % EPSG 3031 WGS84 Antarctic Polar Stereographic
filename = 'IRH_vertical_strain_rates.tif'; % filename code: IRH_LOCATION_TYPEofDATA_TYPEofMODEL_MODELPARAM1_GRIDZISE_TYPEofDATAandFILTER % MODELPARAM1 = basal strain layer thickness (m) / GRIDSIZE = size of gridding resolution in km

% write to file
geotiffwrite(filename,strain_grid_smooth,R,'CoordRefSysCode', coordRefSysCode);

%%