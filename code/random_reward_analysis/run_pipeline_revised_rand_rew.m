clear
clc
close all

%% very similar to the run_pipeline script, with minor changes
% (the laps defined by the IR, not the reward)

frame_rate = 15.49;

belt_length = 170;

data_folder = 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_folder);

file_names = {...
    '9N1_003_000'
    '9N3_003_000'
    '9O2_003_003'...
    };


for i = 1 : length(file_names)
    disp( file_names{i} );
    load([file_names{i}, '_quadrature.mat'], 'quad_data');
    quad_data = double(quad_data);
    frame_num = length(quad_data);
    % save the original frame number (before omitting first frames)
    save([file_names{i}, '_frame_num.mat'], 'frame_num')


    load([file_names{i}, '_DATA.mat'])


    [hardware_ind, quad_data_shifted, calib_value,...
        remove_us_vec, ticks_2_cm] = ...
        get_hardware_data_for_rand_rew(...
        DATA, quad_data, belt_length);

    if isfile([file_names{i}, '_elimination_movie_based.mat'])
        load([file_names{i}, '_elimination_movie_based.mat'],...
            'elimination_movie_based');
        ind_vec = 1 : length(remove_us_vec);
        if isempty(elimination_movie_based.start)
            elimination_movie_based.start = 1;
        end
        if isempty(elimination_movie_based.end)
            elimination_movie_based.end = length(remove_us_vec);
        end
        remove_us_vec(ind_vec < elimination_movie_based.start |...
            ind_vec > elimination_movie_based.end) = true;
    end

    save([file_names{i}, '_hardware_ind.mat'], 'hardware_ind')
    save([file_names{i}, '_calib_value.mat'], 'calib_value')
    save([file_names{i}, '_remove_us_vec.mat'], 'remove_us_vec')

    h_f = plot_trial(frame_num, frame_rate, belt_length,...
        mod(quad_data_shifted, calib_value) * ticks_2_cm,...
        hardware_ind.valve_ind, hardware_ind.lick_ind, hardware_ind.IR_ind,...
        remove_us_vec);
    title(file_names{i}, 'interpreter', 'none')
    drawnow();
    saveas(h_f, [file_names{i}, '_trial.jpeg']);
    saveas(h_f,[file_names{i}, '_trial']);

    start_here = 6;

    % shift everything according to start_here
    quad_data = quad_data(start_here : end);

    field_names = fieldnames(hardware_ind);
    for f = 1 : length(field_names)
        hardware_ind.(field_names{f}) = ...
            hardware_ind.(field_names{f}) - start_here + 1;
        out_of_range = hardware_ind.(field_names{f}) < 1;
        hardware_ind.(field_names{f})(out_of_range) = [];
    end

    quad_data_by_lap = ...
        mod(quad_data - quad_data(hardware_ind.IR_ind(2)), calib_value);
    remove_us_vec = logical(remove_us_vec(start_here : end));


    quad_data_by_lap(remove_us_vec) = nan;

    quad_data_shift =...
        shift_quad_to_fix_IR_bin(quad_data_by_lap, hardware_ind,...
        calib_value);

    if strcmp(file_names{i}(5:7), '003')  % align the track to the familiar rew location
        disp('IR shift by 0.55')
        IR_shift = 0.55;
    end
    quad_bin_number = 10;
    [lap_vec_custom, discrete_quad_data, quad_data_norm] = ...
        create_lap_vec_custom_rand_rew(quad_data, quad_bin_number, hardware_ind,...
        remove_us_vec, IR_shift);

    save([file_names{i}, '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
        'lap_vec_custom', 'discrete_quad_data', 'quad_data_norm',...
        'quad_bin_number');

    save([file_names{i}, '_remove_us_vec_trimmed.mat'],...
        'remove_us_vec');

    %%
    quad_data_smooth = smooth(quad_data, 'moving', 5);
    movement = [0; diff(quad_data_smooth)>1];
    velocity_by_quad_data_smooth = [0; diff(quad_data_smooth)];


    save([file_names{i}, '_quad_data_shift_movement_velocity_hardware.mat'],...
        'quad_data_shift', 'movement',...
        'velocity_by_quad_data_smooth', 'hardware_ind',...
        'start_here');

end



