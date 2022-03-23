%%
clear all
addpath (genpath('D:\R_University_Edinburgh\Toolbox'));
addpath (genpath('D:\R_University_Edinburgh\MacGregor_2016'));

%% Import tif files
tif_1 = 'D:/R_University_Edinburgh/WAIS_accumulation/final_grids/R2_depth_gridded_PIG_IMAFI_THW_smooth_1km_30_18_ALIGNED_IRHsOutline.tif';
tif_2 = 'D:/R_University_Edinburgh/WAIS_accumulation/final_grids/R2_PIG-IMAFI-THW_IceBelow_Fractional_1_GRDSMOOTH_30_18_ALIGNED_IRHsOutline.tif';
tif_3 = 'D:/R_University_Edinburgh/WAIS_accumulation/final_grids/R2_PIG-IMAFI-THW_IceBelow_1_GRDSMOOTH_30_18_ALIGNED_IRHsOutline.tif';

tif_1_vals = geotiffread(tif_1);
tif_2_vals = geotiffread(tif_2);
tif_3_vals = geotiffread(tif_3);

tif_1_vals (tif_1_vals < 0) = NaN;
tif_2_vals (tif_2_vals < 0) = NaN;
tif_3_vals (tif_3_vals < 0) = NaN;

tif_1_vals = flipud(tif_1_vals); 
tif_2_vals = flipud(tif_2_vals); 
tif_3_vals = flipud(tif_3_vals); 

%% extract xy coordinates from grids
tif_info = geotiffinfo(tif_1);

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

%% Import boundaries
drainage = shaperead('D:/R_University_Edinburgh/WAIS_accumulation/calculate_accumulation/IceBoundaries_Antarctica_v2_IRHSOutlines.shp');
xx= []; yy = [];
for i = 1:length(drainage)
    x = drainage(i).X; 
    y = drainage(i).Y;
    xx{i} = x;
    yy {i}= y;
end
xx = horzcat(xx{:}); yy = horzcat(yy{:});
x_ice = xx.'; y_ice = yy.';

IRH_extent = shaperead('D:/R_University_Edinburgh/WAIS_accumulation/calculate_accumulation/R2_WAIS_clipped_extent_final_15032022.shp');
x = IRH_extent.X; y = IRH_extent.Y;
x_IRH = x.'; y_IRH = y.';

%% plot subplot figures
% call figure
f = fullfig; % use fullfig function (Chad Green)

% set plotting variables
plots = {tif_1_vals tif_2_vals tif_3_vals}; % change here if more than 3x plots
axes = ([-1.8e6 -0.6e6 -1e6 5e5]); % set axes

% set titles parameters
titles = {'4.7 ka layer depth (m)' '4.7 ka layer fractional depth' 'Distance from 4.7 ka layer to bed (m)'}; % change here if more than 3x plots
ax = NaN(1, 3); % change here if more than 3x plots
letters = 'a':'c'; % change here if more than 3x plots

% set colormap min/max values
% change here min/max values depending on plot
min_cmap = ([290, 0.15, 270]);
max_cmap = ([2000, 0.8, 2700]);
ytick_cbar_min = ([500, 0.2, 500]);
ytick_cbar_max = ([2000, 0.8, 2500]);

% plot figures
for i = 1:length(plots)
    ax(i) = subplot(1,3,i);
    subaxis(1,3,i,'SpacingVert',0.06,'ML',0.05,'MR',0.03); % change here if more than 3x plots
    axis (axes) % set figure axes
    pcolor(x_grd, y_grd, plots{i}) % raster plot
    shading interp
    caxis([min_cmap(i) max_cmap(i)]);
    axis xy image
    hold on
    plot (x_ice, y_ice,'k','LineWidth',0.3) % ice divides
    hold on
    plot (x_IRH, y_IRH,'k','LineWidth',1) % model boundary
    
    % axes parameters
    if i == 1
        xlabel('PSX (km) ', 'FontSize', 16);
        ylabel('PSY (km) ', 'FontSize', 16);
    end
    axis equal
    xticks([-1.6e6,-1.2e6,-8e5]);
    yticks([-8e5,-6e5,-4e5,-2e5,0e5,2e5]);
    ytickangle(90)
    set(gca,'XTickLabel',[-1600 -1200 -800], 'FontSize', 14);
    set(gca,'YTickLabel',[-800 -600 -400 -200 0 200], 'FontSize', 14);
    
    % colorbar parameters
    cbar = colorbar('fontsize', 14);
    colormap(flipud(viridis));
    set (cbar, 'YTick',linspace(ytick_cbar_min(i),ytick_cbar_max(i),5))   % change here if more than 5x label ticks
    
    % figure parameters
    text(0.0001,1.05,['(' letters(i) ') ' titles{i}],'Units','normalized','FontSize',18);  % change here if more than 3x plots
    grid off
    box on
end

% save figure
% cd 'D:\R_University_Edinburgh\WAIS_accumulation\figures\21032022_figures'
% print('Depth_figure.png','-dpng','-r600');

%%