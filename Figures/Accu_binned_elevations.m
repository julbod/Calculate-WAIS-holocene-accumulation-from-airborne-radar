%% This code calculates accumulation rates binned by elevation
% Note that for the code to run, it is important that all the tif files are
% gridded onto the same grid with the same resolution and size.
%
% Code last updated 15/03/2022 by J. Bodart (UoE)
%
%% Code explanation
% The code imports 3x tif files (1x DEM, 1x modern accumulation, 1x holocene
% accumulation from IRHs) and bins the accumulation values per elevation
% band (e.g. every 50 meters). The code then calculates (a) averaged 
% accumulation rates, (b) total accumulation rates,
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
% Download: https://doi.org/10.5067/E1QL9HFQ7A8M
DEM_tif = 'D:\R_University_Edinburgh\WAIS_accumulation\aligned_grids\BedMachine_surfaceHeight_clipped_1km_final_aligned_IRHsOutline.tif';
[surf, R_DEM, surfs] = geotiffread(DEM_tif);
surf=double(surf);

% load RACMO Surface Mass Balance data product
% Reference: Van Wessem et al., 2014 (RACMO2.3)
% Details: period: 1979-2011; units: m/yr ice equivalent; grid res: 1 km
% Download: https://doi.org/10.1594/PANGAEA.734145     
data_tif = 'D:\R_University_Edinburgh\WAIS_accumulation\aligned_grids\RACMO_accu_clipped_1km_final_aligned_IRHsOutline.tif';
[gridded_RACMO, R_grid, grids] = geotiffread(data_tif);
gridded_RACMO=double(gridded_RACMO);

% load gridded IRH accumulation data from this study
% Reference: Bodart et al., 2022
% Details: 4.72 ka; units: m/yr ice equivalent; grid res: 1 km
% Download: 
IRH_tif = 'D:\R_University_Edinburgh\WAIS_accumulation\aligned_grids\R2_PIG-IMAFI-THW_ACCU_NYE_1_GRDSMOOTH_natural_500m_30_18_ALIGNED_IRHsOutline.tif';
[IRH_vals, R_grid, grids] = geotiffread(IRH_tif);
IRH_vals=double(IRH_vals);

%% plot contours and extract useful information (if required)
% [C,h] = contour(x_grd,y_grd,surf, 0:50:max(surf(:))); % calculate and plot contours
% c_lev=h.LevelList; % get all contour levels
% ind=find(C(1,:)==max(c_lev(c_lev<max(c_lev)))); % find index to the contour level of interest (one below the max in this case)
% n_points=C(2,ind); % number of points along the contour

%% calculate accumulation per elevation bin
% prepare vectors prior to computation
surf_contours = surf(:);
surf_contours(surf_contours<0)=NaN;
gridded_b_RACMO = gridded_RACMO(:);
gridded_b_RACMO(gridded_b_RACMO<0)=NaN;
gridded_bdot = IRH_vals(:);
gridded_bdot(gridded_bdot<0)=NaN;

% round-up surface elevation data to nearest 'cnt' value (if needed)
cnt = 50;
%surf_contours = round(surf_contours/cnt)*cnt; % round to nearest 25 m elevations

% bin all values based on elevation
[Y,E] = discretize(surf_contours,cnt); % Y = bin / E = start value of each bin
YY=discretize(surf_contours,E);
% [c d] = unique(Y); % check if there is any unique elements

% calculate mean accumulation rate for all gridded values within each
% elevation bin
mean_b_RACMO = splitapply(@nanmean,gridded_b_RACMO, YY);
mean_bdot = splitapply(@nanmean,gridded_bdot, YY);

% calculate total accumulation for all gridded values within each
% elevation bin
sum_b_RACMO = splitapply(@nansum,gridded_b_RACMO, YY)./1e4; % convert to Gt
sum_bdot = splitapply(@nansum,gridded_bdot, YY)./1e4; % convert to Gt

% calculate cummulative sum of total accumulation above for all gridded 
% values within each elevation bin
cum_b_RACMO = cumsum(sum_b_RACMO);
cum_bdot = cumsum(sum_bdot);

% calculate x-axis elevation bin
X=[min(surf_contours):(max(surf_contours)-min(surf_contours))/(cnt-1):max(surf_contours)];

%% plot data as 3x vertical plots
figure;
tiledlayout(3,1)

% mean accumulation
ax1 = nexttile;
plot(X,mean_b_RACMO)
hold on
plot(X,mean_bdot)
ylabel('Accum rate (m yr^1)', 'FontSize', 14)
set(gca,'YTickLabel',[0.1 0.2 0.3 0.4 0.5],'FontSize', 14);
set(gca,'XTickLabel',[0 500 1000 1500 2000 2500 3000],'FontSize', 14);
title ('Accumulation rate per elevation', 'FontSize', 18)
legend({'Modern (RACMO)','Holocene (4.72 ka)'},'Location','northeast','Orientation','vertical', 'FontSize', 14)

% total accumulation
ax2 = nexttile;
plot(X,sum_b_RACMO)
hold on
plot(X,sum_bdot)
ylabel('Total accum (Gt yr^1)', 'FontSize', 14)
set(gca,'YTickLabel',[0 0.2 0.4 0.6 0.8 1],'FontSize', 14);
set(gca,'XTickLabel',[0 500 1000 1500 2000 2500 3000],'FontSize', 14);
title ('Total accumulation rate per elevation', 'FontSize', 18)

% cummulative accumulation
ax3 = nexttile;
plot(X,cum_b_RACMO)
hold on
plot(X,cum_bdot)
xlabel('Elevation (m)', 'FontSize', 16)
ylabel('Cummulative accum (Gt yr^1)', 'FontSize', 14)
set(gca,'YTickLabel',[0 5 10 15 20],'FontSize', 14);
set(gca,'XTickLabel',[0 500 1000 1500 2000 2500 3000],'FontSize', 14);
title ('Cummulative accumulation rate per elevation', 'FontSize', 18)
ylim([0 17])

% save plot
% cd 'D:\R_University_Edinburgh\WAIS_accumulation\figures\21032022_figures';
% print(gcf, 'bin_accumulation_graph_v2.png','-dpng','-r600')

%%