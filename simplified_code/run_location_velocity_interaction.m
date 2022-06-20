
pop_mean_STD_loc = nan(size(file_names));
pop_mean_STD_vel = pop_mean_STD_loc;

hist_lims = [-0.5 : 0.01 : 0.5];
hist_vals = nan(1, length(hist_lims) - 1);
vel_bins = [0 : 0.1 : 1];

%% LOAD VARIABLES:
% discrete_location, velocity, movement, event_mat
% cell_IDS - to be included in analysis
%%
velocity(movement~=1) = nan;
velocity(...
    isnan(discrete_location)) = nan;

velocity_discretized = ...
    normalize(velocity', 'range');
velocity_discretized = ...
    discretize(velocity_discretized,...
    vel_bins);


%% samples spent in each condition

num_of_events = ones(size(discrete_location))';
time_spent = ...
    number_of_sim_events_in_bin(discrete_location,...
    velocity_discretized,...
    num_of_events,...
    2);

%% location-velocity heatmaps for single ROIs
loc_vel_interaction_3D_mat = ...
    nan(size(time_spent, 1), size(time_spent, 2), length(cell_IDS));
mean_STD_loc = nan(length(cell_IDS), 1);
mean_STD_vel = mean_STD_loc;

for i = 1 : length(cell_IDS)
    loc_vel_interaction_3D_mat(:, :, i) = ...
        number_of_sim_events_in_bin(discrete_location,...
        velocity_discretized,...
        events_above_min(:, cell_IDS(i)), 1);

    [mean_STD_loc(i), mean_STD_vel(i)] = ...
        STD_loc_vel_interaction(loc_vel_interaction_3D_mat(:, :, i),...
        time_spent);

    % mean_STD_loc: velocity given loc
    % mean_STD_vel: loc given velocity
end


%% location-velocity heatmaps for the population
concurrent_events = sum(events_above_min(:, cell_IDS), 2);
population_mat = ...
    number_of_sim_events_in_bin(discrete_location,...
    velocity_discretized,...
    concurrent_events, 1);

h_im = imagesc(population_mat, 'AlphaData', ~isnan(population_mat));
h_im.Parent.TickDir = 'out';
h_im.Parent.Box = 'off';
colorbar();
xlabel('location');
ylabel('velocity');
[pop_mean_STD_loc, pop_mean_STD_vel] = ...
    STD_loc_vel_interaction(population_mat, time_spent);


%% STD histogram stuff
figure;
h = histogram(mean_STD_vel - mean_STD_loc, hist_lims, ...
    'Normalization', 'probability');
xlabel('mean_STD_vel - mean_STD_loc', 'interpreter', 'none')
ylabel('probability');
ylim([0 0.3])






%% functions:

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
