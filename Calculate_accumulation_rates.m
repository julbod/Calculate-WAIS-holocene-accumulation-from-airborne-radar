%% This code calculates accumulation rates for the 4.72 ka IRH over the WAIS
% see Bodart et al. (2023, The Cryosphere) for more details.
%
% Explanation of code:
% First, the code imports the IRH from a text file, sets-up the model
% parameters, then calculates the accumulation and accumulation uncertainty. 
% The code then grids the values on a 1-km grid, filters the data using a 
% gaussian filter, and exports the data as Geotiff and tabular data files.
% 
% Note that some functions used here originate from the Antarctic Mapping
% Tools toolbox developped by Chad Green (and available at:
% https://github.com/chadagreene/Antarctic-Mapping-Tools)
%
% Code written by J. Bodart (UoE) - 23/02/2022
%
%%
clear all

%% import IRH along-track data from text file
R2 = ['..\IRH_resampled_500m.csv'];
fid = fopen(R2);
R2 = textscan(fid,'%f %f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);

% extract variables from file
R2_x = (R2{1,1}); % PSX
R2_y = (R2{1,2}); % PSY
R2_depth = (R2{1,3}); % Depth of IRH (Fig. 2a)
R2_iceThick = (R2{1,4}); % Ice thickness for each IRH point
R2_strain = (R2{1,5}); % Longitudinal strain rate for each IRH point (Fig. S2b)

% Convert longitudinal strain rate to units equivalent to vertical strain 
% rates (see Fig. S2a)
R2_strain = R2_strain.*1e-4;

% combine individual arrays into one large array
R2_array = horzcat(R2_x,R2_y,R2_depth,R2_iceThick,R2_strain);

%% calculate accumulation rates and uncertainties
for ll = 1:length(R2_array);

    % set up variables for model (ice-core age & age_uncert, thickness, 
    % and depth for each iteration {ll})
    age = vertcat(4.723e3); % Ice-core age for IRH
    age_uncert = vertcat(0.28e3); % Ice-core age-uncertainty for IRH
    thick =  R2_array(ll,4); % Radar ice thickness for IRH
    depth = R2_array(ll,3); % Radar IRH depth
    strain_rate = (R2_array(ll,5)); % Longitudinal strain rate for IRH

    % assign low positive strain value if strain is negative (non-physical)
    if strain_rate < 0
        strain_rate = 1e-6; % low positive value
    end
        
    %% Calculate Nye accumulations (Fig. 3a) and accumulation uncertainties 
    % arising from the age uncertainty of the IRH only, in the first
    % instance
    accum_nye(ll) = deal(-log(1 - (depth ./ thick)) .* (thick ./ age));
    accum_nye_uncert_lower(ll) = deal(-log(1 - (depth ./ thick)) .* (thick ./ (age + age_uncert)));
    accum_nye_uncert_upper(ll) = deal(-log(1 - (depth ./ thick)) .* (thick ./ (age - age_uncert)));

    %% Calculate Shallow-strain accumulations (see Fig. S3) to obtain
    % uncertainty in accumulation rates which are affected by strain rate
    accum_shallow(ll) = (strain_rate*depth*exp(age.*strain_rate))/(exp(age.*strain_rate) - 1);
    accum_shallow_uncert_lower(ll) = (strain_rate*depth*exp((age + age_uncert).*strain_rate))/(exp((age + age_uncert).*strain_rate) - 1);
    accum_shallow_uncert_upper(ll) = (strain_rate*depth*exp((age - age_uncert).*strain_rate))/(exp((age - age_uncert).*strain_rate) - 1);
    
    %% Calculate uncertainty in accumulation for each IRH (Fig S4 a-c)
    lower_uncert_total (ll) = min([accum_nye_uncert_lower(ll) accum_shallow_uncert_lower(ll)]);
    upper_uncert_total (ll) = max([accum_nye_uncert_upper(ll) accum_shallow_uncert_upper(ll)]);
    accum_rel_uncert (ll) = (upper_uncert_total(ll)-lower_uncert_total(ll))/(2.*accum_nye(ll));

    % check progress during loop
    if ll == round(length(R2_array)*0.25)
        display('#### 1/4 of the way there #### ')
    elseif ll == round(length(R2_array)*0.5)
        display('#### 2/4 of the way there #### ')
    elseif ll == round(length(R2_array)*0.75)
        display('#### 3/4 of the way there #### ')
    end
end

%% Export data as tabular array
% concatenate results into array and export as text file (non-gridded data)
IRH_accu = horzcat(R2_array, accum_nye.',lower_uncert_total.',upper_uncert_total.',accum_rel_uncert.');
IRH_accu_combined = table(IRH_accu(:,1),IRH_accu(:,2),IRH_accu(:,3),IRH_accu(:,4),IRH_accu(:,5),IRH_accu(:,6),IRH_accu(:,7),IRH_accu(:,8),IRH_accu(:,9), 'VariableNames', { 'x', 'y','R2_depth','iceThick','longitudinal_strain','Nye_accu','Nye_accu_lower','Nye_accu_upper','Nye_accu_relative'} );
writetable(IRH_accu_combined, '...\accumulation_IRH_combined.txt')

%% Create grid with 1-km spacing around West Antarctica
% centered vertically to get all WAIS IRH coverage
[x,y] = psgrid(-801500, -301500,[2084 1406],1,'xy');

%% Apply filter and grid the along-track data
% smooth along-track variables using moving-average gaussian filter
IRH_accu_smooth = smoothdata(IRH_accu(:,7),'gaussian',30)
IRH_accu_lower_smooth = smoothdata(IRH_accu(:,8),'gaussian',30);
IRH_accu_upper_smooth = smoothdata(IRH_accu(:,9),'gaussian',30);
IRH_accu_rel_uncert_smooth = smoothdata(IRH_accu(:,10),'gaussian',30);

% interpolate filtered along-track data to grid
IRH_accu_grid = griddata(IRH_accu(:,1),IRH_accu(:,2),IRH_accu_smooth,x,y, 'natural');
IRH_accu_lower_grid = griddata(IRH_accu(:,1),IRH_accu(:,2),IRH_accu_lower_smooth,x,y, 'natural');
IRH_accu_upper_grid = griddata(IRH_accu(:,1),IRH_accu(:,2),IRH_accu_upper_smooth,x,y, 'natural');
IRH_accu_rel_uncert_grid = griddata(IRH_accu(:,1),IRH_accu(:,2),IRH_accu_rel_uncert,x,y, 'natural');

% filter grid using average cell filter
filter = fspecial('average',[18,18]);
IRH_accu_grid_smooth = imfilter(IRH_accu_grid,filter);
IRH_accu_lower_grid_smooth = imfilter(IRH_accu_lower_grid,filter);
IRH_accu_upper_grid_smooth = imfilter(IRH_accu_upper_grid,filter);
IRH_accu_rel_uncert_grid_smooth = imfilter(IRH_accu_rel_uncert_grid,filter);

%% export filtered grid data as geotiff file for QGIS import and figures
% get xy limits from grid
xlimits = [min(x(1,:)) max(x(1,:))];
ylimits = [min(y(1,:)) max(y(end,:))];

% specify R object for Geotiff
R_accu = maprefcells(xlimits,ylimits,size(IRH_accu_grid_smooth), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');
R_lower = maprefcells(xlimits,ylimits,size(IRH_accu_lower_grid_smooth), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');
R_upper = maprefcells(xlimits,ylimits,size(IRH_accu_upper_grid_smooth), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');
R_rel = maprefcells(xlimits,ylimits,size(IRH_accu_rel_uncert_grid_smooth), ...
    'ColumnsStartFrom','south','RowsStartFrom','west');

% specify coordinate system and name of files
coordRefSysCode = 3031; % EPSG 3031 WGS84 Antarctic Polar Stereographic
filename_acc = 'IRH_accu.tif';
filename_lower = 'IRH_accu_lower_uncert.tif';
filename_upper = 'IRH_accu_upper_uncert.tif';
filename_rel = 'IRH_accu_relative_uncert.tif';

% write to file
geotiffwrite(filename_acc,IRH_accu_grid_smooth,R_accu,'CoordRefSysCode', coordRefSysCode);
geotiffwrite(filename_lower,IRH_accu_lower_grid_smooth,R_lower,'CoordRefSysCode', coordRefSysCode);
geotiffwrite(filename_upper,IRH_accu_upper_grid_smooth,R_upper,'CoordRefSysCode', coordRefSysCode);
geotiffwrite(filename_rel,IRH_accu_rel_uncert_grid_smooth,R_rel,'CoordRefSysCode', coordRefSysCode);

%%