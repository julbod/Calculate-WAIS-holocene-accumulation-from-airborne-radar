%% This code re-samples IRH depth and ice thickness along 500-m even spacing
% see Bodart et al. (2023, The Cryosphere) for more details.
%
% Explanation of code:
% First, the code imports each IRH file containing the 4.72 ka IRH across
% IMIS, PIG and THW (see paper for details and links for where to access
% these files (all deposited in open-access repositories). The code then 
% finds the distance between each IRH point and a newly created 500-m even 
% spacing line based on the XY information from each IRH vector. It then 
% calcualtes the median IRH depth and ice thickness for all the points that
% fall within 500-m to the right of the point. The code then cleans up any 
% redundant points resulting from calculation of cummulative eucledian 
% distance and exports the data as text file for further processing.
% 
% Note that some functions used here originate from the Antarctic Mapping
% Tools toolbox developped by Chad Green (and available at:
% https://github.com/chadagreene/Antarctic-Mapping-Tools)
%
% Code written by J. Bodart (UoE) - 23/02/2022
%
%%
clear all

%% import IRH data for each individual dataset over IMIS, PIG and THW

% IMIS from Ashmore et al. (2020b)
IMAFI = ['...\Ashmore_IRHs\H2_iceThick_final.txt'];
fid=fopen(IMAFI);
IMAFI=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);
IMAFI_x =(IMAFI{1,1}); % PSX
IMAFI_y =(IMAFI{1,2}); % PSY
IMAFI_depth =(IMAFI{1,3}); % Depth of IRH
IMAFI_iceThick =(IMAFI{1,4}); % Ice thickness at IRH point
IMAFI = horzcat(IMAFI_x,IMAFI_y,IMAFI_depth,IMAFI_iceThick); % combine

% PIG from Bodart et al. (2021b)
PIG =['...\Bodart_IRHs\R2_all_iceThick_final.txt'];
fid=fopen(PIG);
PIG=textscan(fid,'%f %f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);
PIG_x =(PIG{1,1}); % PSX
PIG_y =(PIG{1,2}); % PSY
PIG_depth =(PIG{1,3}); % Depth of IRH
PIG_iceThick =(PIG{1,4}); % Ice thickness at IRH point
PIG = horzcat(PIG_x,PIG_y,PIG_depth,PIG_iceThick); % combine

% THW from Muldoon et al. (2018)
THW =['...\Muldoon_IRHs\LM9_iceThick_final.txt'];
fid=fopen(THW);
THW=textscan(fid,'%f %f %f %f','delimiter',',','headerLines',1);
fclose(fid);
THW_x =(THW{1,1}); % PSX
THW_y =(THW{1,2}); % PSY
THW_depth =(THW{1,3}); % Depth of IRH
THW_depth = THW_depth+10; % add 10 m spatially-invariant firn correction
THW_iceThick =(THW{1,4}); % Ice thickness at IRH point
THW = horzcat(THW_x,THW_y,THW_depth,THW_iceThick); % combine

%% Combine IRH files together
R2_array = vertcat(IMAFI,PIG,THW);

% remove nans in data
nans = all(all(isnan(R2_array),3),2);
R2_array(nans,:,:) = [];

%% Create new 500-m spaced vector based on XY coordinates from IRHs
[xi, yi] = pspath(R2_array(:,1), R2_array(:,2), 500, 'xy'); % re-sample at even points
xy = horzcat(xi,yi); % combine array

%% Re-sample along-track IRH data to 500-m even spacing along flightlines
remove = []; % initialise array
for i = 1:length(xi)
    
    % Compute eucledian distance between the IRH and the 500-m vector 
    % and find values within 500 m of each point
    distance = sqrt(sum(bsxfun(@minus, R2_array(:,1:2), xy(i,1:2)).^2,2));
    bool = distance < 500; % 1 == value within 500-m, 0 == values out
    bool(1:remove) = 0; % remove values 500-m to the left from previous iteration
    
    % Calculate median layer depth and ice thickness along 500-m line
    depth_med (i) = median(R2_array(find(bool==1),3));
    thick_med (i) = median(R2_array(find(bool==1),4));
    
    % Remove used values for next iteration
    remove = length(find(bool == 1));
   
    % Check progress during loop
    if i == round(length(xi)*0.25)
        display('#### 1/4 of the way there #### ')
    elseif i == round(length(xi)*0.5)
        display('#### 2/4 of the way there #### ')
    elseif i == round(length(xi)*0.75)
        display('#### 3/4 of the way there #### ')
    end
end

%% Combine arrays and remove nans arising from pspath calculations
R2_resampled = horzcat(xi,yi,depth_med.', thick_med.');
nans = isnan(R2_resampled(:,3));
R2_resampled(nans==1,:) = []; % remove nans

% export to tabular file
% table = table(R2_resampled(:,1),R2_resampled(:,2),R2_resampled(:,3),R2_resampled(:,4), 'VariableNames', { 'PSX', 'PSY','IRH_depth','iceThick'} );
% writetable(table, '...\R2_combined_resampled_500m.txt')

%% Remove points which are not exactly positioned on original flight line
for i = 1:length(R2_resampled)
    % if closest distance between two sets of XY is larger than safe value
    % (e.g. safe value assumed to be ~ 20-m for max uncertainty in along-track)
    if min(sqrt(sum(bsxfun(@minus, R2_array(:,1:2), R2_resampled(i,1:2)).^2,2))) > 20
        idx (i) = nan; % place nans
    end
    
    % Check progress during loop
    if i == round(length(R2_resampled)*0.25)
        display('#### 1/4 of the way there #### ')
    elseif i == round(length(R2_resampled)*0.5)
        display('#### 2/4 of the way there #### ')
    elseif i == round(length(R2_resampled)*0.75)
        display('#### 3/4 of the way there #### ')
    end
end

% Remove rows of data where minimum distance between the two datasets is
% larger than value specified in loop (e.g. safe value assumed to be 20 m)
R2_resampled(isnan(idx.'),:) = [];

%% Export to tabular file
table = table(R2_resampled(:,1),R2_resampled(:,2),R2_resampled(:,3),R2_resampled(:,4), 'VariableNames', { 'PSX', 'PSY','IRH_depth','iceThick'} );
writetable(table, '...\IRH_resampled_500m.txt')

%%   