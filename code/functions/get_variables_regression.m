function [df_F, events_above_min, movement,...
    quad_data_norm, lap_vec_custom, ...
    velocity_by_quad_data_smooth] = ...
    get_variables_regression(file_name)

load([file_name, '_full_sig_ROIs_df_F.signals'], '-mat', 'df_F');
load([file_name, '_new_event_det_df_F.mat'], 'events_above_min');

load([file_name, '_remove_us_vec_trimmed.mat'], 'remove_us_vec');
load([file_name, '_remove_ROIs.mat'], 'remove_ROIs');
if isfile([file_name, '_elimination_movie_based.mat'])
    load([file_name, '_elimination_movie_based.mat'],...
        'elimination_movie_based');
    remove_ROIs = union(remove_ROIs, elimination_movie_based.dead_ROIs);
end

load([file_name, '_quad_data_shift_movement_velocity_hardware.mat'],...
    'movement', 'hardware_ind',...
    'velocity_by_quad_data_smooth');
load([file_name, '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
    'quad_data_norm', 'lap_vec_custom');

% removing partial laps:
remove_us_vec(isnan(quad_data_norm)) = 1;


df_F(remove_us_vec, :) = [];
events_above_min(remove_us_vec, :) = [];
df_F(:, remove_ROIs) = [];
events_above_min(:, remove_ROIs) = [];

hardware_ind.valve_ind = ...
    setdiff(hardware_ind.valve_ind, find(remove_us_vec==1));
hardware_ind.lick_ind = ...
    setdiff(hardware_ind.lick_ind, find(remove_us_vec==1));
hardware_ind.IR_ind	 = ...
    setdiff(hardware_ind.IR_ind, find(remove_us_vec==1));

movement(remove_us_vec) = [];
velocity_by_quad_data_smooth(remove_us_vec) =[];
quad_data_norm(remove_us_vec) = [];
lap_vec_custom(remove_us_vec) = [];


end
