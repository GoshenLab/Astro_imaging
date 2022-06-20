function [quad_data_shift, shift_by, h_f] = ...
    shift_quad_to_fix_IR_bin(quad_data_by_lap, hardware_ind,...
    calib_value)

% quad_data_by_lap needs to be with all values. NANs according to the
% remove us vector

rew_tick_num_lap = quad_data_by_lap(hardware_ind.IR_ind);
h_f = figure;
subplot(2,2,1)
plot(rew_tick_num_lap, '.')
ylim([0 calib_value])
subplot(2,2,3)
polarscatter(rew_tick_num_lap./calib_value * 2 * pi,...
    ones(size(rew_tick_num_lap)), '.')

% shift everything to be in the zero region
angle_ticks = rew_tick_num_lap./calib_value * 2 * pi;
angle_ticks = repmat(angle_ticks, numel(angle_ticks),1);
angle_diff = mod((angle_ticks - angle_ticks'), 2 * pi);
angle_diff = min(angle_diff, 2 * pi - angle_diff);
shift_by = max(angle_diff(:))/(2 * pi) * calib_value;
quad_data_shift = mod((shift_by + quad_data_by_lap), calib_value);
rew_tick_num_lap = quad_data_shift(hardware_ind.IR_ind);

% make sure the minimal is at zero (there's a safety margin using the
% previous calculation)
quad_data_shift = mod(quad_data_shift - min(rew_tick_num_lap), calib_value);
rew_tick_num_lap = quad_data_shift(hardware_ind.IR_ind);

subplot(2,2,2)
plot(rew_tick_num_lap, '.')
ylim([0 calib_value])
subplot(2,2,4)
polarscatter(rew_tick_num_lap./calib_value * 2 * pi,...
    ones(size(rew_tick_num_lap)), '.')
end
