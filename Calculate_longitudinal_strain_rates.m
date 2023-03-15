%% This code was originally written by Joe MacGregor (NASA Goddard) and
% adapted by Julien Bodart (UoE) to calculate the longitudinal strain rates
% (or velocity gradient/variation) based on present surface velocity, see 
% Bodart et al. (2023, The Cryosphere) for more details.
%
% For details on the original method, refer to Waddington et al. 
% (2007; JGlac); MacGregor et al. (2009; Ann. Glac.) and MacGregor et al. 
% (2016; Science)
%
% Explanation of code:
% First, the code imports modern velocities from MeASURE dataset, resampled
% onto a common 1-km grid. The code then filters the data, fills holes and 
% recalculates speeds where these are less than 100 m a-1, calculates
% longitudinal strain rates, and finally grids and exports the data as 
% geotiff files for import into QGIS and to plot as figures. This data is
% shown in Fig. S2b of Bodart et al. (2023).
%
% Code originally written by J. MacGregor (NASA)
% Code adapted by J. Bodart (UoE) - 23/02/2022
%
%%
clear all

%% Load reference grid (MEASURES V2 ice-flow resampled at a 1-km grid)
tif = '...\MEaSUREs_iceFlow_velocities_1km.tif';
tif_info = geotiffinfo(tif);

% Extract PSX and PSY coordinates from tif file
[x_min, y_min] = deal(tif_info.BoundingBox(1, 1), tif_info.BoundingBox(1, 2));
[num_x, num_y] = deal(tif_info.Width, tif_info.Height);
tif_inc = tif_info.GeoTIFFTags.ModelPixelScaleTag(1);

% Assign coordinate data
x_coords = (x_min + (tif_inc / 2)) + (0:tif_inc:((num_x - 1) * tif_inc));
y_coords = (y_min + (tif_inc / 2)) + (0:tif_inc:((num_y - 1) * tif_inc))';

% Grid PSX and PSY coordinates
[x_grd, y_grd] = meshgrid(x_coords, y_coords);
[num_y, num_x] = size(x_grd);

%% Set parameters for filtering
age_ref                     = 4.72e3; % reference age for layer in years
decim                       = 5; % decimation indices but also effectively km
filt_type                   = 'exp'; % exponential filter
thick_filt_ref              = 10; % number of ice thicknesses over which to run filter

%% Load filtered velocity data from D calculations
load ..\D_results.mat vars_filt

% Decimate regular grids into new grid variables (x_filt, y_filt)
length_grd                  = 1e3 * decim; % grid size, m
ind_x                       = find(~mod(round(x_grd(1, :)), decim), 1):decim:find(~mod(round(x_grd(1, :)), decim), 1, 'last');
ind_y                       = find(~mod(round(y_grd(:, 1)), decim), 1):decim:find(~mod(round(y_grd(:, 1)), decim), 1, 'last');
[x_filt, y_filt]            = deal(x_grd(ind_y, ind_x)./1e3, y_grd(ind_y, ind_x)./1e3);

% Get grid size
[num_decim_y, num_decim_x]  = deal(length(ind_y), length(ind_x));

% Unpack arrays
speed_x_filt = double(vars_filt{2}); speed_y_filt = double(vars_filt{3});
speed_x_uncert_filt = double(vars_filt{6}); speed_y_uncert_filt = double(vars_filt{7});
surf_filt = double(vars_filt{1});

% Make sure data is of type double and not single
x_filt = double(x_filt); y_filt = double(y_filt);

%% Recalculate flow direction for slower ice-flow regions (~ <100 m/yr)
% Calculate magnitude of velocity and its uncertainty
speed_filt                  = sqrt((speed_x_filt .^ 2) + (speed_y_filt .^ 2)); % velocity magnitude prior to filtering
speed_uncert_filt           = sqrt((speed_x_uncert_filt .^ 2) + (speed_y_uncert_filt .^ 2)); % velocity uncertainty magnitude

% Calculate surface-velocity and elevation-gradient flow azimuths
az_speed                    = atan2(speed_y_filt, speed_x_filt); % azimuth of surface velocity
[elev_grad_x, elev_grad_y]  = gradient(surf_filt, length_grd, length_grd); % compute approximate surface elevation gradients
az_elev                     = atan2(-elev_grad_y, -elev_grad_x); % azimuth of gradient in surface elevation

% Smooth azimuths, extract useful elements
[az_sin_cat, az_cos_cat]    = deal(NaN(num_decim_y, num_decim_x, 2));
az_sin_cat(:, :, 1)         = sin(az_speed);
az_cos_cat(:, :, 1)         = cos(az_speed);
az_sin_cat(:, :, 2)         = sin(az_elev);
az_cos_cat(:, :, 2)         = cos(az_elev);

% Weight filter exponentially using reference speed
speed_az_decay              = 100; % speed above which weight is unity for InSAR surface speeds,  m/s
speed_uncert_rel_decay      = 0.1;
wt_az_elev                  = exp(-speed_filt ./ speed_az_decay) + exp(-speed_uncert_rel_decay ./ (speed_uncert_filt ./ speed_filt)); % characteristic length of surface speed to unity ratio
wt_az_elev(wt_az_elev > 1)  = 1;
wt_az_speed                 = 1 - wt_az_elev;
wt_az_elev(isnan(az_speed)) = 1; % maximum weight (1) if the other is NaN
wt_az_speed(isnan(az_elev)) = 1;
az_mean                     = atan2(((az_sin_cat(:, :, 1) .* wt_az_speed) + (az_sin_cat(:, :, 2) .* wt_az_elev)), ((az_cos_cat(:, :, 1) .* wt_az_speed) + (az_cos_cat(:, :, 2) .* wt_az_elev)));
az_sin                      = sin(az_mean);
az_cos                      = cos(az_mean);
az_sin_sign                 = sign(az_sin);
az_cos_sign                 = sign(az_cos);
az_sin_abs_rel              = abs(az_sin) ./ (abs(az_sin) + abs(az_cos));
az_cos_abs_rel              = abs(az_cos) ./ (abs(az_sin) + abs(az_cos));

% Reproject speeds
[speed_x_filt, speed_y_filt]= deal((speed_filt .* az_cos), (speed_filt .* az_sin));
speed_filt                  = sqrt((speed_x_filt .^ 2) + (speed_y_filt .^ 2)); % assign newly filtered velocity magnitude values post filtering
speed_filt_ref              = speed_filt;

%% Compute longitudinal strain rates
[speed_grad_x, speed_grad_y] = gradient(speed_filt); % gradient of filtered surface speed (1/yr)
[speed_x_norm, speed_y_norm] = deal((speed_x_filt ./ speed_filt), (speed_y_filt ./ speed_filt)); % unit vectors for x/y speeds (dimensionless)
speed_grad_lon = (speed_grad_x .* speed_x_norm) + (speed_grad_y .* speed_y_norm); % du/dx (gradient * unit vector, 1/yr)

%% export filtered grid data as geotiff file for QGIS import and figures
% get xy limits from grid
xlimits = [min(x_filt(1,1:end))*1e3 max(x_filt(1,1:end))*1e3]; % in m
ylimits = [min(y_filt(1:end,1))*1e3 max(y_filt(1:end,1))*1e3]; % in m

% specify R object for Geotiff
R = maprefcells(double(xlimits),double(ylimits),size(speed_grad_lon), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');

% specify coordinate system and name of files
coordRefSysCode = 3031; % EPSG
filename = 'IRH_longitudinal_strain_rates.tif'; 

% write to file
geotiffwrite(filename,speed_grad_lon,R,'CoordRefSysCode', coordRefSysCode);

%%