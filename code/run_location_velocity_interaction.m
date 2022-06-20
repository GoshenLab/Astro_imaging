
clear
clc
% close all


plot_loc_tuning_curves = false;
plot_vel_tuning_curves = false;
plot_location_velocity_ROI_heatmaps = false;


file_names = {...
    '8C3_100_001',...
    '9B2_100_004',...
    '9H3_003_000',...
    '9I3_004_000',...
    '9P2_000_000',...
    '9N1_000_000',...
    '9Q4_000_000',...
    '9O3_000_000',...
    '9S3_000_000',...
    };

data_folder = 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_folder);

pop_mean_STD_loc = nan(size(file_names));
pop_mean_STD_vel = pop_mean_STD_loc;

hist_lims = [-0.5 : 0.01 : 0.5];
hist_vals = nan(length(file_names), length(hist_lims) - 1);
vel_bins = [0 : 0.1 : 1];

for f = 1 : length(file_names)

    file_name = file_names{f};
    
    load([file_name '_new_event_det_df_F.mat'], 'events_above_min');
    
    cell_IDS = 1 : size(events_above_min, 2);
    
    load([file_name, '_remove_ROIs.mat'], 'remove_ROIs');
    if isfile([file_name, '_elimination_movie_based.mat'])
        load([file_name, '_elimination_movie_based.mat'],...
            'elimination_movie_based');
        remove_ROIs = union(remove_ROIs, elimination_movie_based.dead_ROIs);
    end
    
    cell_IDS = setdiff(cell_IDS, remove_ROIs)';
        
    
    load([file_name, '_remove_us_vec_trimmed.mat'], 'remove_us_vec');
    
    
    load([file_name, '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
        'discrete_quad_data');
   
    load([file_name, '_quad_data_shift_movement_velocity_hardware.mat'],...
        'velocity_by_quad_data_smooth', 'movement')
    
    
    velocity_by_quad_data_smooth(movement~=1) = nan;
    velocity_by_quad_data_smooth(...
        isnan(discrete_quad_data)) = nan; % to include the remove_us_vec
    
    velocity_by_quad_data_smooth_discretized = ...
        normalize(velocity_by_quad_data_smooth', 'range');
    velocity_by_quad_data_smooth_discretized = ...
        discretize(velocity_by_quad_data_smooth_discretized,...
        vel_bins);
    
    
    %% samples spent in each condition
    
    num_of_events = ones(size(discrete_quad_data))';
    
    time_spent = ...
        number_of_sim_events_in_bin(discrete_quad_data,...
        velocity_by_quad_data_smooth_discretized,...
        num_of_events,...
        2);
    
       
    %% location-velocity heatmaps for single ROIs
    loc_vel_interaction_3D_mat = ...
        nan(size(time_spent, 1), size(time_spent, 2), length(cell_IDS));
    mean_STD_loc = nan(length(cell_IDS), 1);
    mean_STD_vel = mean_STD_loc;
    
    for i = 1 : length(cell_IDS)
        loc_vel_interaction_3D_mat(:, :, i) = ...
            number_of_sim_events_in_bin(discrete_quad_data,...
            velocity_by_quad_data_smooth_discretized,...
            events_above_min(:, cell_IDS(i)), 1);
        
        [mean_STD_loc(i), mean_STD_vel(i)] = ...
            STD_loc_vel_interaction(loc_vel_interaction_3D_mat(:, :, i),...
            time_spent);
        
        % mean_STD_loc: velocity given loc
        % mean_STD_vel: loc given velocity
    end


    %% population mat
    concurrent_events = sum(events_above_min(:, cell_IDS), 2);
    population_mat = ...
        number_of_sim_events_in_bin(discrete_quad_data,...
        velocity_by_quad_data_smooth_discretized,...
        concurrent_events, 1);
    
    h_im = imagesc(population_mat, 'AlphaData', ~isnan(population_mat));
    h_im.Parent.TickDir = 'out';
    h_im.Parent.Box = 'off';
    colorbar();
    xlabel('location');
    ylabel('velocity');
    [pop_mean_STD_loc(f), pop_mean_STD_vel(f)] = ...
        STD_loc_vel_interaction(population_mat, time_spent);
    
    title({file_name;...
        'mean concurrent events as a function of velocity and location';...
        'for the entire population';...
        ['mean STD_loc = ', num2str(pop_mean_STD_loc(f)), '     ', ...
        'mean STD_vel = ', num2str(pop_mean_STD_vel(f))]},...
        'interpreter', 'none')
    
    %% STD histogram stuff
    figure;
    h = histogram(mean_STD_vel - mean_STD_loc, hist_lims, ...
        'Normalization', 'probability');
    hist_vals(f, :) = h.Values;
    title(file_name, 'interpreter', 'none')
    xlabel('mean_STD_vel - mean_STD_loc', 'interpreter', 'none')
    ylabel('probability');
    ylim([0 0.3])
end

%%
figure;

subplot(1,3,1)
plot(hist_lims(1 : end-1) + 0.005, hist_vals', 'LineWidth', 2);
line([0 0],[0 0.3], 'linestyle' , '--','color', 'k')
legend(file_names, 'Interpreter', 'none')
xlabel('mean_STD_vel - mean_STD_loc', 'interpreter', 'none')
ylabel('probability');

subplot(1,3,2)
shadedErrorBar(hist_lims(1 : end-1) + 0.005, mean(hist_vals),...
    std(hist_vals)./sqrt(length(file_names)));
line([0 0],[0 0.3], 'linestyle' , '--','color', 'k')
xlabel('mean_STD_vel - mean_STD_loc', 'interpreter', 'none')
ylabel('probability');


%%
h_f = figure;
plot(ones(size(pop_mean_STD_loc)), pop_mean_STD_loc, '.', 'MarkerSize', 20);
hold on;
plot(ones(size(pop_mean_STD_vel))*2, pop_mean_STD_vel, '.', 'MarkerSize', 20);
plot([pop_mean_STD_loc; pop_mean_STD_vel], 'Color', [.5 .5 .5])
xlim([0.5 2.5])
ylim([0 0.13])
h_f.Children.TickDir = 'out';
h_f.Children.Box = 'off';
ylabel('Mean Weighted STD')


%% mean STD of row/column of interaction matrix
% time_spent(isnan(time_spent)) = 0;

function [mean_STD_loc, mean_STD_vel] = ...
    STD_loc_vel_interaction(data_mat, time_spent)
STD_loc = nan(size(time_spent, 2), 1);
STD_vel = nan(size(time_spent, 1), 1);

for i = 1 : size(time_spent, 2)
    STD_loc(i) = std(data_mat(:, i), time_spent(:, i), 'omitnan');
end

for i = 1 : size(time_spent, 1)
    STD_vel(i) = std(data_mat(i, :), time_spent(i, :), 'omitnan');
end

mean_STD_loc = mean(STD_loc);
mean_STD_vel = mean(STD_vel);
end
%% functions:

function plot_lap_and_pooled_tuning_curves(cell_IDS, sig_3D_mat,...
    pooled_mean_activity, x_label, y_label)

figure;
c = 1;
ROWS = 6;
COLS = 10;

for i = 1 : length(cell_IDS)
    curr_data = sig_3D_mat(:, :, cell_IDS(i))';
    
    subplot(ROWS, COLS, c)
    imagesc(curr_data, 'AlphaData', ~isnan(curr_data));
    xlabel(x_label)
    ylabel(y_label)
    title(['ROI #', num2str(cell_IDS(i))]);
    
    subplot(ROWS, COLS, c + COLS)
    shadedErrorBar([], mean(curr_data, 'omitnan'),...
        std(curr_data, 'omitnan')./sqrt(sum(~isnan(curr_data))),...
        'lineprops', 'b');
    hold on;
    plot(pooled_mean_activity(:, cell_IDS(i)), 'r');
    xlabel(x_label)
    ylabel('mean activity probability')
    if c == ROWS * COLS - COLS
        figure;
        c = 1;
    else
        if mod(c, COLS) == 0
            c = c + 1 + COLS;
        else
            c = c + 1;
        end
    end
    
end
end