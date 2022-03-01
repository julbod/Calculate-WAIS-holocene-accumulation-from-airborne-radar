%% This code calculates accumulation rates binned by elevation
% Note that for the code to run, it is important that all the tif files are
% gridded onto the same grid with the same resolution and size. Ideally,
% this can be done in QGIS using the 'Raster align' tool. 
%
% Code last updated 01/03/2022 by J. Bodart (UoE)
%
%% Code explanation
% The code imports 3x tif files (1x DEM, 1x modern accumulation, 1x holocene
% accumulation from IRHs), then rounds-up the elevation values to nearest
% desired value, and calculates the elevation bins. The code then
% calculates (a) averaged accumulation rates, (b) total accumulation rates,
% and (c) cummulative accumulation rates, that are found within each
% elevation bin. Finally, the code plots the data and saves it as png file.
%
%% set variables and directories
clear all
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));

%% load grids 
% load REMA surface height data from BedMachine data product (modified)
% Reference: Howat et al., 2019 (REMA); Morlighem et al., 2019 (BedMachine)
% Details: nominal date: 2012; units: m (height); grid res: 1 km
% Download: https://doi.org/10.5067/E1QL9HFQ7A8MDEM_tif = 'D:\R_University_Edinburgh\WAIS_accumulation\aligned_grids\BedMachine_surfaceHeight_clipped_1km_final_aligned_IRHsOutline.tif';
[surf, R_DEM, surfs] = geotiffread(DEM_tif);
surf (surf==0) = NaN;
surf=double(surf);

% load RACMO snow accumulation data from ALBMAP data product
% Reference: VanDeBerg et al., 2006 (RACMO) & Le Brocq et al., 2010 (ALBMAP) 
% Details: period: 1980-2004; units: m/yr ice equivalent; grid res: 1 km
% Download: https://doi.org/10.1594/PANGAEA.734145  
data_tif = 'D:\R_University_Edinburgh\WAIS_accumulation\aligned_grids\ALBMAP_accu_clipped_1km_final_aligned_IRHsOutline.tif';
[gridded_vals, R_grid, grids] = geotiffread(data_tif);
gridded_vals (gridded_vals==0) = NaN;
gridded_vals=double(gridded_vals);

% load gridded IRH accumulation data from this study
% Reference: Bodart et al., 2022
% Details: 4.72 ka; units: m/yr ice equivalent; grid res: 1 km
% Download: 
IRH_tif = 'D:\R_University_Edinburgh\WAIS_accumulation\aligned_grids\R2_PIG-IMAFI-THW_ACCU_NYE_1_GRDSMOOTH_ALIGNED_IRHsOutline.tif';
[IRH_vals, R_grid, grids] = geotiffread(IRH_tif);
IRH_vals (IRH_vals==0) = NaN;
IRH_vals=double(IRH_vals);

%% plot contours and extract useful information (if required)
% [C,h] = contour(x_grd,y_grd,surf, 0:50:max(surf(:))); % calculate and plot contours
% c_lev=h.LevelList; % get all contour levels
% ind=find(C(1,:)==max(c_lev(c_lev<max(c_lev)))); % find index to the contour level of interest (one below the max in this case)
% n_points=C(2,ind); % number of points along the contour

%% calculate accumulation per elevation bin
% prepare vectors prior to computation
surf_contours = surf(:);
gridded_b = gridded_vals(:);
gridded_bdot = IRH_vals(:);

% round-up surface elevation data to nearest 'cnt' value
cnt = 25;
surf_contours = round(surf_contours/cnt)*cnt; % round to nearest 25 m elevations

% bin all values based on elevation
[Y,E] = discretize(surf_contours,cnt); % Y = bin / E = start value of each bin
YY=discretize(surf_contours,E);
% [c d] = unique(Y); % check if there is any unique elements

% calculate mean accumulation rate for all gridded values within each
% elevation bin
mean_b = splitapply(@nanmean,gridded_b, YY);
mean_bdot = splitapply(@nanmean,gridded_bdot, YY);

% calculate total accumulation for all gridded values within each
% elevation bin
sum_b = splitapply(@nansum,gridded_b, YY)./1e4; % convert to Gt
sum_bdot = splitapply(@nansum,gridded_bdot, YY)./1e4; % convert to Gt

% calculate cummulative sum of total accumulation above for all gridded 
% values within each elevation bin
cum_b = cumsum(sum_b);
cum_bdot = cumsum(sum_bdot);

% calculate x-axis elevation bin
X=[min(surf_contours):(max(surf_contours)-min(surf_contours))/(cnt-1):max(surf_contours)];

%% plot data as 3x vertical plots
tiledlayout(3,1)

% mean accumulation
ax1 = nexttile;
plot(X,mean_b)
hold on
plot(X,mean_bdot)
ylabel('Accumulation rate (m yr^1)')
title ('Accumulation rate per elevation')
legend({'Modern (RACMO)','Holocene (4.72 ka)'},'Location','northeast','Orientation','vertical')

% total accumulation
ax2 = nexttile;
plot(X,sum_b)
hold on
plot(X,sum_bdot)
ylabel('Total accumulation (Gt yr^1)')
title ('Total accumulation rate per elevation')

% cummulative accumulation
ax3 = nexttile;
plot(X,cum_b)
hold on
plot(X,cum_bdot)
xlabel('Elevation (m)')
ylabel('Cummulative accumulation (Gt yr^1)')
title ('Cummulative accumulation rate per elevation')
ylim([0 17])

% save plot
% cd 'D:\R_University_Edinburgh\WAIS_accumulation\figures';
% print(gcf, 'bin_accumulation_graph.png','-dpng','-r600')

%%