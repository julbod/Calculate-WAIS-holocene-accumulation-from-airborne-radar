%% This code was adapted from Joe MacGregor's (NASA) code to calculate the
% suitability of the Local-Layer Approximation for calculating accumulation
% rates in englacial layers using 1-D modelling.

% For details on the method, refer to Waddington et al. (2007; JGlac);
% MacGregor et al. (2009; Ann. Glac.) and MacGregor et al. (2016; Science)

% Code last updated 17/02/2021 by J. Bodart (UoE)

%% Code explanation
% The gridded boundary condition inputs should all be of similar resolution
% If grid resolution is more than half an order of magnitude difference,
% then users should interpolate the grid to reconcile grid cell resolution.
% Ideally, all grids should have the same resolution and be mapped
% onto the same grid (1 km gridding resolution, as used here, is ideal).
%
% Coordinate system used here is WGS84: EPSG 3031 and units are in m.
%
% Note that the REMA surface elevation product was downloaded from
% BedMachine NetCDF, and converted from surface height (in ice equivalent)
% to the original REMA surface product following the equation: 
% Z = surf + firn + geoid.
%
% Once the different data products are inputed, 
% the code filters the data using a thickness-dependent exponential 
% decaying filter. Depending on size of study area, it is advised to only
% do this once and save the output variables as .mat files to save on 
% processing power and time. If the code throws an error towards the end 
% of the computation, it may be due to grid resolution and lengths of 
% the input products (if those where not re-gridded previously).
%
% The code then recalculates flow direction for slower ice flow regions
% (e.g. interior; < ~100 m/yr) where surface velocity directions are
% noisier than the signal due to ionospheric noise which limits the quality
% of satellite-derived surface velocity there (see MacGregor et al., 2016;
% Science, S.I.). The code combines the azimuths inferred from the modern
% ice-flow velocities with the gradients in present surface elevations
% to infer the surface-velocity field there based on both the observed
% surface velocity and its relative uncertainty.
%
% The code finally uses pool computing (not necessary however can speed things
% up if well done) to calculate the Lpath, Lb, Lh and finally the D number.
% This can then be plotted geographically to identify areas where the LLA
% may not be suitable and 1-D model results should not be trusted.

%% set variables and directories
clear all
cd 'D:\R_University_Edinburgh\MEaSUREs_ice_velocity';
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));

%% load coordinates from reference grid (i.e. MEASURES V2 ice-flow 1-km grid)
tif = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\MEASURES_speed_clipped_1km_final.tif';
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

%% set parameters for filtering
age_ref                     = 4.72e3; % reference age for layer in kyr
decim                       = 5; % decimation, indices but also effectively km
filt_type                   = 'exp'; % triang or exp
scale_type                  = 'thick'; % ind or thick
ind_filt                    = 5; % width of filter, indices, used only when scale_type=ind, should be an odd integer
thick_filt_ref              = 10; % number of ice thicknesses over which to run filter, used only when scale_type=thick

%% load data products and filter them using thickness-dependent filter
do_filter = false; % default is not to run following code and use saved .mat file
if do_filter % however, if 'do_filter = true', run loop
    disp('Importing boundary condition datasets...')
    
    %% load ice-flow velocities from InSAR MEaSUREs V2 data product 
    % Reference: Rignot et al., 2017
    % Details: period: 1996-2016; units: m/yr; grid res: 1 km
    % Download: https://doi.org/10.5067/D7GK8F5J8M8R
    tif_VX = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\MEASURES_speedVX_clipped_1km_final.tif'; % X velocity components
    tif_VY = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\MEASURES_speedVY_clipped_1km_final.tif'; % Y velocity components
    tif_ERRX = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\MEASURES_speedUncertainty_clipped_1km_final.tif'; % Error X velocity components
    tif_ERRY = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\MEASURES_speedUncertaintyERRY_clipped_1km_final.tif'; % Error Y velocity components
    
    % extract data from GeoTifs
    speed_x = geotiffread(tif_VX);
    speed_y = geotiffread(tif_VY);
    speed_x_uncert = geotiffread(tif_ERRX);
    speed_y_uncert = geotiffread(tif_ERRY);
    
    % remove NaNs and unrealistic negative values
    speed_x(speed_x < -3620) = NaN;
    speed_y(speed_y < -3850) = NaN;
    speed_x_uncert(speed_x_uncert < 0) = NaN;
    speed_y_uncert(speed_y_uncert < -0.4) = NaN;
    
    % flip array in up/down direction
    speed_x = flipud(speed_x);
    speed_y = flipud(speed_y);
    speed_x_uncert = flipud(speed_x_uncert);
    speed_y_uncert = flipud(speed_y_uncert);
    
    % calculate velocity magnitudes
    speed_xy = sqrt((speed_x .^ 2) + (speed_y .^ 2));
    speed_xy (speed_xy > 4200) = NaN; % remove large outliers (numbers based on tiff in QGIS)
    speed_xy_uncert = sqrt((speed_x_uncert .^ 2) + (speed_y_uncert .^ 2));
    
    %% load RACMO snow accumulation data from ALBMAP data product
    % Reference: VanDeBerg et al., 2006 (RACMO) & Le Brocq et al., 2010 (ALBMAP) 
    % Details: period: 1980-2004; units: m/yr ice equivalent; grid res: 1 km
    % Download: https://doi.org/10.1594/PANGAEA.734145    
    tif = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\ALBMAP_accu_clipped_1km_final.tif';
    precip = geotiffread(tif);
    precip(precip < -10) = NaN; % remove large values
    precip(precip < 0) = 0; % replace small outliers smaller than 0 to 0
    
    % flip array in up/down direction
    precip = flipud(precip);
    
    %% load ice thickness data from BedMachine data product
    % Reference: Morlighem et al., 2019 
    % Details: nominal date: 2012; units: m; grid res: 1 km
    % Download: https://doi.org/10.5067/E1QL9HFQ7A8M
    tif = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\BedMachine_thick_clipped_1km_final.tif';
    thick = geotiffread(tif);
    
    % flip array in up/down direction
    thick = flipud(thick);
    
    %% load REMA surface height data from BedMachine data product (modified)
    % Reference: Howat et al., 2019 (REMA); Morlighem et al., 2019 (BedMachine)
    % Details: nominal date: 2012; units: m (height); grid res: 1 km
    % Download: https://doi.org/10.5067/E1QL9HFQ7A8M
    tif = 'D:\University_Edinburgh\QGIS_Linux\Chapt_3_accumulation\clipped_rasters\BedMachine_surfaceHeight_clipped_1km_final.tif';
    surf = geotiffread(tif);
    
    % flip array in up/down direction
    surf = flipud(surf);
    
    %% plot data (if necessary)
    % decim = 10; % decimate matrices to speed up plotting
    % var_x = x_coords; var_y = y_coords; var_val = speed_x; % specify vars
    % figure('position', [5 5 750 750])
    % hold on
    % pcolor(var_x(1:decim:end, 1:decim:end), var_y(1:decim:end, 1:decim:end), double(var_val(1:decim:end, 1:decim:end)))
    % shading interp
    % caxis([min(var_val(:)) max(var_val(:))])
    % axis xy image
    % xlabel('Polar stereographic X (km)')
    % ylabel('Polar stereographic Y (km)')
    % colorbar('fontsize', 20)
    % grid on
    % box on

    %% apply filter to data
    disp('Starting filtering of datasets...')
    
    % set up filtering grid array based on ice thickness
    filt_decim                  = cell(1, max([1 round(1e-3 * max(thick(:)) * thick_filt_ref)]));
    for ii = 1:length(filt_decim)
        filt_decim{ii}          = repmat(exp([-(ii:-1:1) 0 -(1:ii)] ./ thick_filt_ref), ((2 * ii) + 1), 1) .* repmat(exp([-(ii:-1:1) 0 -(1:ii)] ./ thick_filt_ref)', 1, ((2 * ii) + 1));
        filt_decim{ii}          = filt_decim{ii} ./ sum(filt_decim{ii}(:));
    end
    
    % give name to non-filtered and filtered variables and place in giant
    % cell arrays (vars = non-filted; vars_filt = fully filtered)
    vars                        = {'surf' 'speed_x' 'speed_y' 'thick' 'precip' 'speed_x_uncert' 'speed_y_uncert'};
    vars_filt                   = {'surf_filt' 'speed_x_filt' 'speed_y_filt' 'thick_filt' 'precip_filt' 'speed_x_uncert_filt' 'speed_y_uncert_filt'};
    num_var                     = length(vars);
    
    % decimate regular grids into new grid variables (x_filt, y_filt)
    length_grd                  = 1e3 * decim; % grid size, m
    ind_x                       = find(~mod(x_grd(1, :), decim), 1):decim:find(~mod(x_grd(1, :), decim), 1, 'last');
    ind_y                       = find(~mod(y_grd(:, 1), decim), 1):decim:find(~mod(y_grd(:, 1), decim), 1, 'last');
    [x_filt, y_filt]            = deal(x_grd(ind_y, ind_x)./1e3, y_grd(ind_y, ind_x)./1e3); % XY grid for new interpolated grid
    
    % get grid size
    [num_decim_y, num_decim_x]  = deal(length(ind_y), length(ind_x));
    
    % set up giant cell array by filling it with nans
    for ii = 1:num_var
        vars_filt{ii} = NaN(num_decim_y, num_decim_x);
    end
    
    % place all non-filtered gridded variables into giant cell array
    vars = cell(1,numel(vars_filt));
    vars{1} = surf;
    vars{2} = speed_x;
    vars{3} = speed_y;
    vars{4} = thick;
    vars{5} = precip;
    vars{6} = speed_x_uncert;
    vars{7} = speed_y_uncert;
    
    % reference set of ice thicknesses (equivalent to 10 H)
    thick_ref = (1:length(filt_decim)) ./ (1e-3 * thick_filt_ref);
    
    % find closest matching set of ice thicknesses
    filt_decim_ref = interp1(thick_ref, 1:length(filt_decim), thick(ind_y, ind_x), 'nearest', 'extrap');
    
    % start filtering all gridded products and re-grid to same grid
    for jj = 1:num_decim_x
        disp([num2str(jj) '/' num2str(num_decim_x)])
        if ~mod(jj, 50)
            disp(jj)
        end
        
        for ii = 1:num_decim_y
            % current filter matrix and half-size
            filt_decim_curr     = filt_decim{filt_decim_ref(ii, jj)};
            num_filt_curr       = (size(filt_decim_curr, 1) - 1) / 2;
            
            % indices within original matrix
            ind_tmp             = -num_filt_curr:num_filt_curr;
            ii_tmp              = repmat((ind_y(ii) + ind_tmp)', ((2 * num_filt_curr) + 1), 1);
            jj_tmp              = ind_x(jj) + ind_tmp(ones(((2 * num_filt_curr) + 1), 1), :);
            jj_tmp              = jj_tmp(:);
            ind_good            = find((ii_tmp > 0) & (jj_tmp > 0) & (ii_tmp <= num_y) & (jj_tmp <= num_x));
            ind_curr            = sub2ind([num_y num_x], ii_tmp(ind_good), jj_tmp(ind_good));
            
            % filter each variable
            for kk = 1:num_var
                vars_filt{kk}(ii, jj) = nansum(filt_decim_curr(ind_good) .* vars{kk}(ind_curr)) * ((numel(filt_decim_curr) ^ 2) / (length(ind_good) * length(find(~isnan(vars{kk}(ind_curr))))));
            end
        end
    end
    
    % save cell array to .mat file for convenience
    disp('Saving newly filtered data to .mat file...')
    cd 'D:\R_University_Edinburgh\WAIS_accumulation\calculate_D_parameter';
    save('filter_vars_v2.mat','vars_filt')
    
% else if 'do_filter = false'
else
    disp('Loading already filtered datasets...')
    % load .mat variable containing all filtered and re-gridded data
    load D:\R_University_Edinburgh\WAIS_accumulation\calculate_D_parameter\filter_vars_v2.mat vars_filt
    
    % decimate regular grids into new grid variables (x_filt, y_filt)
    length_grd                  = 1e3 * decim; % grid size, m
    ind_x                       = find(~mod(x_grd(1, :), decim), 1):decim:find(~mod(x_grd(1, :), decim), 1, 'last');
    ind_y                       = find(~mod(y_grd(:, 1), decim), 1):decim:find(~mod(y_grd(:, 1), decim), 1, 'last');
    [x_filt, y_filt]            = deal(x_grd(ind_y, ind_x)./1e3, y_grd(ind_y, ind_x)./1e3);
    
    % get grid size
    [num_decim_y, num_decim_x]  = deal(length(ind_y), length(ind_x));
end

%% unpack arrays
surf_filt = double(vars_filt{1}); speed_x_filt = double(vars_filt{2}); speed_y_filt = double(vars_filt{3});
thick_filt = double(vars_filt{4}); precip_filt = double(vars_filt{5});
speed_x_uncert_filt = double(vars_filt{6}); speed_y_uncert_filt = double(vars_filt{7});

% make sure data is of type double and not single
x_filt = double(x_filt); y_filt = double(y_filt);

%% Recalculate flow direction for slower ice-flow regions (~ <100 m/yr)
% calculate magnitude of velocity and its uncertainty
speed_filt                  = sqrt((speed_x_filt .^ 2) + (speed_y_filt .^ 2)); % velocity magnitude prior to filtering
speed_uncert_filt           = sqrt((speed_x_uncert_filt .^ 2) + (speed_y_uncert_filt .^ 2)); % velocity uncertainty magnitude

% calculate surface-velocity and elevation-gradient flow azimuths
az_speed                    = atan2(speed_y_filt, speed_x_filt); % azimuth of surface velocity
[elev_grad_x, elev_grad_y]  = gradient(surf_filt, length_grd, length_grd); % compute approximate surface elevation gradients
az_elev                     = atan2(-elev_grad_y, -elev_grad_x); % azimuth of gradient in surface elevation

% smooth azimuths, extract useful elements
[az_sin_cat, az_cos_cat]    = deal(NaN(num_decim_y, num_decim_x, 2));
az_sin_cat(:, :, 1)         = sin(az_speed);
az_cos_cat(:, :, 1)         = cos(az_speed);
az_sin_cat(:, :, 2)         = sin(az_elev);
az_cos_cat(:, :, 2)         = cos(az_elev);

% weight filter exponentially using reference speed
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

% reproject speeds
[speed_x_filt, speed_y_filt]= deal((speed_filt .* az_cos), (speed_filt .* az_sin));
speed_filt                  = sqrt((speed_x_filt .^ 2) + (speed_y_filt .^ 2)); % assign newly filtered velocity magnitude values post filtering
speed_filt_ref              = speed_filt;

%% load or calculate D parameter
do_D = false; % default is not to run following code and use saved .mat file
if do_D % however, if 'do_filter = true', run loop

    %% calculate D parameter
    disp('Starting D calculation...')

    % initializations
    [num_decim_y, num_decim_x] = size(speed_filt);
    [length_path, length_precip, length_thick] = deal(NaN(num_decim_y, num_decim_x));
    warning('off', 'MATLAB:polyfit:PolyNotUnique')

    for ii = 1:num_decim_y
        disp([num2str(ii) '/' num2str(num_decim_y)])

        % parfor jj = 1:num_decim_x % use pool computing
        for jj = 1:num_decim_x

            % if any nans, unable to compute values and thus skip iteration
            if any(isnan([speed_filt(ii, jj) thick_filt(ii, jj) precip_filt(ii, jj)]))
                continue
            end

            % compute reversed flowlines using stream2 function
            % "A streamline can be thought of as the path a massless particle 
            % takes flowing through a velocity field"
            fl = stream2(x_filt, y_filt, -speed_x_filt, -speed_y_filt, x_filt(ii, jj), y_filt(ii, jj), [0.5 1000]); % reversed flowline
            fl = fl{1};

            % if any nans, trim them
            if any(isnan(fl(:, 1)))
                fl = fl(1:(find(isnan(fl(:, 1)), 1) - 1), :);
            end

            % if flowline went nowhere after this, skip iteration
            if (size(fl, 1) < 2)
                continue
            end

            % assign flowline XY values
            [x_fl, y_fl] = deal(fl(:, 1)', fl(:, 2)');

            % linearly interpolate speed along flowline
            speed_fl = interp2(x_filt, y_filt, speed_filt_ref, x_fl, y_fl, '*linear', NaN);

            % fix nans in speed data if any
            if any(isnan(speed_fl))
                [x_fl, y_fl, speed_fl] = deal(x_fl(1:(find(isnan(speed_fl), 1) - 1)), y_fl(1:(find(isnan(speed_fl), 1) - 1)), speed_fl(1:(find(isnan(speed_fl), 1) - 1)));
            end

            % if flowline went nowhere after this, skip iteration
            if (length(x_fl) < 2)
                continue
            end

            % calculate distance and age of particule along flowline
            dist_fl = cumsum([0 sqrt((diff(x_fl) .^ 2) + (diff(y_fl) .^ 2))]); % distance along flowline
            age_fl = cumsum([0 (diff(1e3 .* dist_fl) ./ speed_fl(1:(end - 1)))]); % age along flowline

            % calculate horizontal distance of particule from surface to
            % position of current layer in ice (a.k.a. "Lpath" parameter)
            length_path(ii, jj) = interp1(age_fl, dist_fl, age_ref, 'linear', 'extrap'); % length of path to reach current layer

            % remove any values which are larger than the reference age of the
            % layer of interest
            if any(age_fl > age_ref)
                [x_fl, y_fl, dist_fl] = deal(x_fl(1:(find((age_fl > age_ref), 1) - 1)), y_fl(1:(find((age_fl > age_ref), 1) - 1)), dist_fl(1:(find((age_fl > age_ref), 1) - 1))); % trim flowline coordinates to reference age
            end

            % if flowline went nowhere after this, skip iteration
            if (length(x_fl) < 2)
                continue
            end

            % end flowline at Lpath
            [x_fl, y_fl, dist_fl] = deal([x_fl interp1(dist_fl, x_fl, length_path(ii, jj), 'linear', 'extrap')], [y_fl interp1(dist_fl, y_fl, length_path(ii, jj), 'linear', 'extrap')], [dist_fl length_path(ii, jj)]);

            % calculate "Lbdot" and "H" parameters
            precip_fl = interp2(x_filt, y_filt, precip_filt, x_fl, y_fl, '*linear', NaN); % Lbdot
            thick_fl = interp2(x_filt, y_filt, thick_filt, x_fl, y_fl, '*linear', NaN); % H

            % calculate "dbdot/dx" and "dH/dx"
            poly_precip = polyfit(dist_fl, precip_fl, 1); % linear regression, y = (dbdot/dx)*x + b
            poly_thick = polyfit(dist_fl, thick_fl, 1); % y = (dH/dx)*x + b

            % calculate characteristic lengths "Lb" and "Lh" parameters
            length_precip(ii, jj) = abs(poly_precip(1) ./ nanmean(precip_fl)); % 1/L_bdot = |(1/bdot)*dbdot/dx|
            length_thick(ii, jj) = abs(poly_thick(1) ./ nanmean(thick_fl)); % 1/L_H = |(1/H)*dH/dx|
        end
    end
    warning('on', 'MATLAB:polyfit:PolyNotUnique')

    % calculate D number and re-assign data to final arrays
    % (unitless: <1 == LLA likely valid; >1 == LLA likely invalid)
    D = length_path .* (length_precip + length_thick);
    [length_precip, length_thick] = deal((1 ./ length_precip), (1 ./ length_thick));

    %% save and export data
    disp('Done calculating D, now saving...')

    cd ('D:\R_University_Edinburgh\WAIS_accumulation\calculate_D_parameter');
    save(['D_' num2str(1e-3 * age_ref) 'ka.mat'], '-v7.3', 'length_path', 'length_precip', 'length_thick', 'D', 'num_decim_y', 'num_decim_x', 'age_ref')

    disp(['Saved file'])
    %delete(pool)

% else if 'do_D = false'
% extract all data points equal or smaller than threshold in D array
else
    disp('Loading already calculated D...')
    cd ('D:\R_University_Edinburgh\WAIS_accumulation\calculate_D_parameter');
    load 'D_4.72ka.mat'   
    
    %% mask data
    D_mask = D;
    threshold = 1.25; % set threshold: <1 == LLA likely valid; >1 == LLA likely invalid
    idx = find(D_mask > threshold); % indices higher than threshold
    D_mask(idx) = NaN; % replace by nans
    
    %% contour data
    xy_D05                      = contourc(x_filt(1, :), y_filt(:, 1), D, [0.5 0.5]);
    [x_D05, y_D05]              = deal(xy_D05(1, :), xy_D05(2, :));
    %save mat/xy_D05 x_D05 y_D05

    xy_D1                       = contourc(x_filt(1, :), y_filt(:, 1), D, [1 1]);
    [x_D1, y_D1]                = deal(xy_D1(1, :), xy_D1(2, :));
    %save mat/xy_D1 x_D1 y_D1
    
    % figure; contour(x_filt(1, :), y_filt(:, 1), D, 2);
end

%% write data to NetCDF
% cd 'D:\R_University_Edinburgh\WAIS_accumulation\D_code_breakdown\test\post-filter'; % specify directory
% var_id = 4; % specify which variable you want in large array

% % define a netcdf
% ncid = netcdf.create('BedMachine_surface_filtereddddd.nc','NOCLOBBER'); % first create netcdf
% dim_x = netcdf.defDim(ncid, 'X-coordinates', length(x_filt(1,1:end))); % specify amount of dimensions
% dim_y = netcdf.defDim(ncid, 'Y-coordinates', length(y_filt(1:end,1))); % specify amount of dimensions
% varid_1 = netcdf.defVar(ncid,'X-coordinates','NC_DOUBLE',dim_x); % specify type of byte
% varid_2 = netcdf.defVar(ncid,'Y-coordinates','NC_DOUBLE',dim_y); % specify type of byte
% varid_3 = netcdf.defVar(ncid,'surface','NC_DOUBLE',[dim_x,dim_y]); % specify type of byte
% netcdf.endDef(ncid); % come out of define mode
%
% % write to netcdf
% netcdf.putVar(ncid,varid_1,x_filt(1,1:end));
% netcdf.putVar(ncid,varid_2,y_filt(1:end,1));
% netcdf.putVar(ncid,varid_3,vars_filt{var_id});
% netcdf.close(ncid) % close netcdf
%
% % read netcdf to check
% % ncdisp('filtered_regridded_thicknesses.nc')
%
%% write data to GeoTIF
% get xy limits from gridded map
% xlimits = [min(x_filt(1,1:end))*1e3 max(x_filt(1,1:end))*1e3]; % in m
% ylimits = [min(y_filt(1:end,1))*1e3 max(y_filt(1:end,1))*1e3]; % in m
% 
% % specify R object
% R = maprefcells(double(xlimits),double(ylimits),size(D), ...
%     'ColumnsStartFrom','south','RowsStartFrom','west');
% 
% coordRefSysCode = 3031; % EPSG
% filename = 'D_4.72ka.tif'; 
% geotiffwrite(filename,D,R,'CoordRefSysCode', coordRefSysCode);

%%