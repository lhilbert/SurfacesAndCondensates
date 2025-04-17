clc
clear all
global conditions_folders conditions_fields num_images num_ROIs ROI_size centered_area_size

tic  

% Specify folders
base_folders_with_surface = { 
    '', 
    ''
};

base_folders_without_surface = {
    '',
    ''
};

num_images = 10; % Number of images per condition 
num_ROIs = 10; % Number of desired ROIs
ROI_size = 100; % ROI size
centered_area_size = 700; % Central zone size

% Initialization
conditions_folders = {'0_0ul','0_5ul', '0_75ul', '1_0ul', '1_25ul', '1_5ul', '1_75ul', '2_0ul'};
conditions_fields = {'ul_0_0', 'ul_0_5', 'ul_0_75', 'ul_1_0', 'ul_1_25', 'ul_1_5', 'ul_1_75', 'ul_2_0'};
results_with_surface = struct();
results_without_surface = struct();

for i = 1:length(conditions_fields)
    results_with_surface.(conditions_fields{i}) = struct('mean_condensate_intensity', [], ...
                                                       'mean_background_intensity', [], ...
                                                       'num_condensates', [], ...
                                                       'number_of_valid_ROI', 0);  
    
    results_without_surface.(conditions_fields{i}) = struct('mean_condensate_intensity', [], ...
                                                          'mean_background_intensity', [], ...
                                                          'num_condensates', [], ...
                                                          'number_of_valid_ROI', 0);  
end

% Processing of data with surface
for base_idx = 1:length(base_folders_with_surface)
    results_with_surface = process_datasets(base_folders_with_surface{base_idx}, results_with_surface);
end
% Save
save(fullfile(base_folders_with_surface{1}, 'results_with_surface_full.mat'), 'results_with_surface');

disp('Process done - with surface.');

% Processing of data without surface
for base_idx = 1:length(base_folders_without_surface)
    results_without_surface = process_datasets(base_folders_without_surface{base_idx}, results_without_surface);
end
% Save
save(fullfile(base_folders_without_surface{1}, 'results_without_surface_full.mat'), 'results_without_surface');

disp('Process done - without surface.');
%% 

% Plot
plot_results_bootstrap(results_with_surface, results_without_surface);
%% 

elapsed_time = toc;  
fprintf('Elapsed time : %.2f secondes\n', elapsed_time);

function results = process_datasets(base_folder, results)
    global conditions_folders conditions_fields num_images num_ROIs ROI_size centered_area_size
    min_size = 4; % 'Cleaning' of 2x2 objects
    
    for cond_idx = 1:length(conditions_folders)
        folder_name = conditions_folders{cond_idx};
        field_name = conditions_fields{cond_idx};
        input_folder = fullfile(base_folder, folder_name);
        
        if ~exist(input_folder, 'dir')
            continue;  % Skip if folder doesn't exist
        end
        
        roi_folder = fullfile(input_folder, 'ROI_centralregion');
        
        if exist(roi_folder, 'dir')
            delete(fullfile(roi_folder, '*.png'));
        else
            mkdir(roi_folder);
        end

        number_of_valid_ROI = 0; 
        
        for i = 1:num_images
            raw_image_path = fullfile(input_folder, 'segmentation_stdDev4_full_tif', sprintf('zstackSlice_Condensate_raw_%d.tif', i));
            mask_path = fullfile(input_folder, 'segmentation_stdDev4_full', sprintf('zstackSlice_Condensate_seg_mask_%d.png', i));
            
            if ~isfile(raw_image_path) || ~isfile(mask_path)
                continue;
            end
            
            fprintf('\nImage %d in %s\n', i, input_folder);

            raw_image = double(imread(raw_image_path)); 
            mask = imread(mask_path);
            mask = bwareaopen(mask, min_size);
            mask = mask > 0; 

            se = strel('disk', 2);  
            mask = imclose(mask, se); 

            mask = imfill(mask, 'holes');
            
            mask = uint8(mask) * 255;
            
            [height, width] = size(raw_image);
            
            % Define central zone
            center_x = floor(width/2);
            center_y = floor(height/2);
            x_min = max(1, center_x - centered_area_size/2);
            x_max = min(width, center_x + centered_area_size/2);
            y_min = max(1, center_y - centered_area_size/2);
            y_max = min(height, center_y + centered_area_size/2);
            
            for j = 1:num_ROIs
                
                x = randi([x_min, x_max - ROI_size + 1]);
                y = randi([y_min, y_max - ROI_size + 1]);
                
                roi_raw = raw_image(y:y+ROI_size-1, x:x+ROI_size-1);
                roi_mask = mask(y:y+ROI_size-1, x:x+ROI_size-1);
                
                imwrite(uint16(roi_raw), fullfile(roi_folder, sprintf('ROI_raw_%d_%d.png', i, j)));
                imwrite(roi_mask, fullfile(roi_folder, sprintf('ROI_mask_%d_%d.png', i, j)));
                
                % Mean condensate intensity
                condensate_pixels = roi_raw(roi_mask > 0);
                if isempty(condensate_pixels)
                    mean_condensate_intensity = NaN;
                else
                    mean_condensate_intensity = mean(condensate_pixels(:));
                    number_of_valid_ROI = number_of_valid_ROI + 1; 
                end
                
                background_pixels = roi_raw(roi_mask == 0);
                mean_background_intensity = mean(background_pixels(:));

                CC = bwconncomp(roi_mask);
                num_condensates = CC.NumObjects;

                % Stock results
                results.(field_name).mean_condensate_intensity = [results.(field_name).mean_condensate_intensity; mean_condensate_intensity];
                results.(field_name).mean_background_intensity = [results.(field_name).mean_background_intensity; mean_background_intensity];
                results.(field_name).num_condensates = [results.(field_name).num_condensates; num_condensates];
            end
        end

                results.(field_name).number_of_valid_ROI = number_of_valid_ROI;
    end
end

function plot_results_bootstrap(results_with_surface, results_without_surface)
    conditions = fieldnames(results_with_surface);
    concentrations = {'0','0.18','0.27','0.36','0.45','0.54','0.63','0.72'};
    figure('Position', [100, 100, 1600, 400]);
    
    means_with = zeros(1, length(conditions));
    means_without = zeros(1, length(conditions));
    errors_with = zeros(1, length(conditions));
    errors_without = zeros(1, length(conditions));
    num_condensates_with = zeros(1, length(conditions));
    num_condensates_without = zeros(1, length(conditions));
    
    bg_means_with = zeros(1, length(conditions));
    bg_means_without = zeros(1, length(conditions));
    bg_errors_with = zeros(1, length(conditions));
    bg_errors_without = zeros(1, length(conditions));

    for i = 1:length(conditions)
        field = conditions{i};
        
        data_with = results_with_surface.(field).mean_condensate_intensity;
        data_without = results_without_surface.(field).mean_condensate_intensity;
        
        data_with = data_with(~isnan(data_with) & ~isinf(data_with));
        data_without = data_without(~isnan(data_without) & ~isinf(data_without));

        bg_with = results_with_surface.(field).mean_background_intensity;
        bg_without = results_without_surface.(field).mean_background_intensity;

        bg_with = bg_with(~isnan(bg_with) & ~isinf(bg_with));
        bg_without = bg_without(~isnan(bg_without) & ~isinf(bg_without));
        
        num_condensates_with(i) = length(data_with);
        num_condensates_without(i) = length(data_without);

        total_condensates_with(i) = sum(results_with_surface.(field).num_condensates, 'omitnan');
        total_condensates_without(i) = sum(results_without_surface.(field).num_condensates, 'omitnan');
        number_of_ROI_with(i) = results_with_surface.(field).number_of_valid_ROI;
        number_of_ROI_without(i) = results_without_surface.(field).number_of_valid_ROI;
        
        means_with(i) = mean(data_with, 'omitnan');
        means_without(i) = mean(data_without, 'omitnan');

        if num_condensates_with(i) > 1
            errors_with(i) = std(data_with) / sqrt(length(data_with));
            disp(std(data_with))
        end
        if num_condensates_without(i) > 1
            errors_without(i) = std(data_without) / sqrt(length(data_without));
            disp(std(data_without))
        end

        bg_means_with(i) = mean(bg_with, 'omitnan');
        bg_means_without(i) = mean(bg_without, 'omitnan');
        
        if ~isempty(bg_with) && length(bg_with) > 1
            bg_errors_with(i) = std(bg_with) / sqrt(length(bg_with));
            disp(std(bg_with)) 
        end
        
        if ~isempty(bg_without) && length(bg_without) > 1
            bg_errors_without(i) = std(bg_without) / sqrt(length(bg_without));
            disp(std(bg_without)) 
        end

    end
    
    % Marker size
    base_size = 25;  % Base size
    number_of_ROI_with(number_of_ROI_with == 0) = 1; 
    number_of_ROI_without(number_of_ROI_without == 0) = 1;
    number_of_stacks = 200;

    marker_sizes_with = base_size * total_condensates_with./number_of_stacks;
    marker_sizes_without = base_size * total_condensates_without./number_of_stacks;

    marker_sizes_with = max(marker_sizes_with, 1);
    marker_sizes_without = max(marker_sizes_without, 1);
    
    % Plot with surface
    subplot(2, 2, 1);
    errorbar([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_with, errors_with, 'b-', 'LineWidth', 1.5, 'capsize', 10);
    hold on;
    errorbar([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], bg_means_with, bg_errors_with, bg_errors_with, 'ko-', 'LineWidth', 1.5, 'capsize', 10);
    scatter([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_with, marker_sizes_with, 'b');
    xlabel('Concentration (µM)'); ylabel('Intensity');
    title('With Surface'); grid on;
    % set(gca, 'XTick', [0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], 'XTickLabel', concentrations);
    legend('Condensates','Background', 'Location', 'best');
    
    % Plot without surface
    subplot(2, 2, 3);
    errorbar([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_without, errors_without, 'r-', 'LineWidth', 1.5, 'capsize', 10);
    hold on;
    errorbar([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], bg_means_without, bg_errors_without, bg_errors_without, 'ko-', 'LineWidth', 1.5, 'capsize', 10);
    scatter([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_without, marker_sizes_without, 'r', 'filled');
    xlabel('Concentration (µM)'); ylabel('Intensity');
    title('Without Surface'); grid on;
    % set(gca, 'XTick', [0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], 'XTickLabel', concentrations);
    legend('Condensates','Background', 'Location', 'best');
    
    % Compare
    subplot(2, 2, [2 4]);
    errorbar([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_with, errors_with, 'bo-', 'LineWidth', 1.5);
    hold on;
    errorbar([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_without, errors_without, 'ro-', 'LineWidth', 1.5);
    scatter([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_with, marker_sizes_with, 'b', 'filled');
    scatter([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], means_without, marker_sizes_without, 'r', 'filled');
    xlabel('Concentration (µM)'); ylabel('Mean Condensate Intensity');
    title('Comparison'); grid on;
    % set(gca, 'XTick', [0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], 'XTickLabel', concentrations);
    legend('With Surface', 'Without Surface', 'Location', 'best');
    
    % % Plot - number of condensates
    % subplot(1, 4, 4);
    % bar([0, 0.18, 0.27, 0.36, 0.45, 0.54, 0.63, 0.72], [total_condensates_with; total_condensates_without]');
    % xlabel('Concentration (µM)'); ylabel('Number of Condensates (total count)');
    % title('Condensate Count'); grid on;
    % set(gca, 'XTick', 1:length(conditions), 'XTickLabel', concentrations);
    % legend('With Surface', 'Without Surface', 'Location', 'best');
    
end
