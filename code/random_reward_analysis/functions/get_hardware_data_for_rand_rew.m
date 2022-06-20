function [hardware_ind, quad_data_shifted, calib_value, remove_us_vec, ticks_2_cm] = ...
    get_hardware_data_for_rand_rew(DATA, quad_data, belt_length)

% I updated this on 24.7.20 such that when I remove the bad laps, I still
% leave the relevant reward locations, to enable the custom lap division
% for random reward, the IR ind replaces the valve...

sample_num = 100; % this works when a lap takes longer than 100/sampling rate, i.e. 6.5 sec
valve_ind = ...
    find_events_with_large_ISI(DATA.valve, 0.5, sample_num);
IR_ind = ...
    find_events_with_large_ISI(DATA.IR, 0.5, sample_num);

sample_num_lick = 0;
lick_ind =...
    find_events_with_large_ISI(DATA.lickometer, 0.5, sample_num_lick);

hardware_ind.valve_ind = valve_ind;
hardware_ind.IR_ind = IR_ind;
hardware_ind.lick_ind = lick_ind;

%% get a sense of how uniform the laps are...

IR_tick_diff = diff(quad_data(IR_ind));
figure;
plot(IR_tick_diff, '.')
hold on
plot(1 : length(IR_tick_diff),...
    mean(IR_tick_diff) * ones(1, length(IR_tick_diff)))

%% which laps should be left?
% assuming you get laps that are too short/long - don't use them for
% calculation of the calibration value.

THRESHOLD = 400;
IR_log_vec = ...
    find(abs(mean(IR_tick_diff) - IR_tick_diff) > THRESHOLD);
calib_value = ...
    mean(IR_tick_diff(~(abs(mean(IR_tick_diff) - IR_tick_diff) > ...
    THRESHOLD)));
ticks_2_cm = belt_length / calib_value;

remove_us_vec = zeros(size(quad_data));

for i = 1 : length(IR_log_vec)
    % if it's the last lap:
    if IR_log_vec(i) == length(IR_ind) - 1
% this way you include the valve ind, just not what's after it        
curr_ind = IR_ind(IR_log_vec(i)) + 1 : length(quad_data);
    % if it's the first lap:
    elseif IR_log_vec(i) == 1 
% this way you include the valve ind, just not what's before it        
         curr_ind = 1 : IR_ind(IR_log_vec(i) + 1) - 1;
    % for all other cases:
    else
        curr_ind = ...
            [IR_ind(IR_log_vec(i)) + 1 : IR_ind(IR_log_vec(i) + 1)];
% leave the edge, so you know where the reward was given
    end
    remove_us_vec(curr_ind) = 1;
end

%% create similar quad_data locations across files (just for ploting!!!)
SHIFT_TICKS = 900;
% I chose the second valve as my reference because the first is often
% wrong...
first_lap = quad_data(valve_ind(2));
quad_data_shifted = quad_data - first_lap - SHIFT_TICKS;

end


function index = find_events_with_large_ISI(data, threshold_sig, min_ISI)

a = data > threshold_sig;
[~,index] = findpeaks(double(a), 'MinPeakDistance', min_ISI);

end