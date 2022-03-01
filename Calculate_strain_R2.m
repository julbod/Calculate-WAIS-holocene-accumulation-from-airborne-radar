%% This code calculates vertical strain rates from Nye-calculated accumulation
% rates for R2/H2/LM9 over the WAIS
% 
% First, the code imports the text file containing the accumulation rates
% calculated using 'R2_Nye_accu_gridding.m' file for R2/H2/LM9 over the
% WAIS, and calculates vertical strain rates using simple equation. The
% code finally grids and smoothes the data using gaussian filtering and
% exports the data as TIF and text files.
% 
% Code written by J. Bodart (UoE) - 28/02/2022
%
%%
clear all
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));

%% import layer data
R2_accu = ['D:\R_University_Edinburgh\WAIS_accumulation\calculate_accumulation\Line_accumulation\R2_PIG-IMAFI-THW_ACCU_NYE_LINE.txt'];
fid=fopen(R2_accu);
R2_accu=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);

x_accu =(R2_accu{1,1});
y_accu =(R2_accu{1,2});
z_depth =(R2_accu{1,3});
iceThick =(R2_accu{1,4});
accu =(R2_accu{1,5});

%% Calculate vertical strain stress
strain = (-(accu))./iceThick;

%% create grid with 1-km spacing around West Antarctica
%[x,y] = psgrid(-77.34,-95.56,600,1,'xy'); % Centered over PIG - resolution grid: 1km
%[x,y] = psgrid(-77.5820,-108.4458,[250 420],5,'xy'); % Centered vertically to get PIG and THW Bodart and Muldoon IRHs - resolution is 5 km grid
%[x,y] = psgrid(-79.5820,-93.56,[1300 750],1,'xy'); % Centered horizontally to get PIG and IMAFI Bodart and Ashmore IRHs - resolution is 1 km grid
%[x,y] = psgrid(-77.3807,-112.4346,[650 750],5,'xy'); % Centered around Thwaites for Leuro layers - resolution is 1 km grid
[x,y] = psgrid(-801500, -301500,[2084 1406],1,'xy'); % Centered vertically to get Bodart, Ashmore and Muldoon IRHs - resolution is 1 km grid

%% create interpolation from scatter with original XY data
% fit interpolated data to grid
strain_grid = griddata(x_accu,y_accu,strain,x,y);

% smooth data using gaussian filter
strain_grid_smooth = smoothdata(strain_grid,'gaussian',20);

%% export grid as geotiff file
cd 'D:\R_University_Edinburgh\WAIS_accumulation\calculate_strain\Gridded_strain'

% get xy limits from gridded map
xlimits = [min(x(1,:)) max(x(1,:))];
ylimits = [min(y(1,:)) max(y(end,:))];

% specify R object
R = maprefcells(xlimits,ylimits,size(strain_grid_smooth), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');

coordRefSysCode = 3031; % EPSG
filename = 'R2_PIG-IMAFI-THW_STRAIN_1_GRDSMOOTH.tif'; % filename code: IRH_LOCATION_TYPEofDATA_TYPEofMODEL_MODELPARAM1_GRIDZISE_TYPEofDATAandFILTER % MODELPARAM1 = basal strain layer thickness (m) / GRIDSIZE = size of gridding resolution in km
geotiffwrite(filename,strain_grid_smooth,R,'CoordRefSysCode', coordRefSysCode);

%% export line data to csv
table = table(x_accu, y_accu, z_depth, iceThick, accu, strain, 'VariableNames', {'x','y','R2_depth','iceThick','Nye_accu','strain'} );
writetable(table, 'D:\R_University_Edinburgh\WAIS_accumulation\calculate_strain\Line_strain\R2_PIG-IMAFI-THW_STRAIN_LINE.txt')

%%