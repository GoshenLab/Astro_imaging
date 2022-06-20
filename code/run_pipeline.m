
frame_rate = 15.49;
belt_length = 170;

% insert data_dir! this code does not work for the random case!
data_dir = 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_dir)
all_file_names = dir('*_quadrature.mat');
all_quad_file_name = {all_file_names.name};
all_file_names = cellfun(@(x) x(1:11), all_quad_file_name, 'uniformoutput', false);
mouse_names = cellfun(@(x) x(1:3), all_file_names, 'uniformoutput', false);

for i = 1 : length(all_file_names)
    file_name = all_file_names{i};
    if sum(strcmp(file_name, {'9N1_003_000', '9N3_003_000', '9O2_003_003'})>0) % don't use this for the random reward!
        continue
    end
    load(all_quad_file_name{i});
    quad_data = double(quad_data);
    frame_num = length(quad_data);

    % save the original frame number (before omitting first frames)
    save([file_name, '_frame_num.mat'], 'frame_num')
    load([file_name, '_DATA.mat'], 'DATA')

    [hardware_ind, quad_data_shifted, calib_value,...
        remove_us_vec, ticks_2_cm] = ...
        get_hardware_data(DATA, quad_data, belt_length);

    if isfile([file_name, '_elimination_movie_based.mat'])
        load([file_name, '_elimination_movie_based.mat'],...
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

    save([file_name, '_hardware_ind.mat'], 'hardware_ind')
    save([file_name, '_calib_value.mat'], 'calib_value')
    save([file_name, '_remove_us_vec.mat'], 'remove_us_vec')

    h_f = plot_trial(frame_num, frame_rate, belt_length,...
        mod(quad_data_shifted, calib_value) * ticks_2_cm,...
        hardware_ind.valve_ind, hardware_ind.lick_ind, hardware_ind.IR_ind,...
        remove_us_vec);
    title(file_name, 'interpreter', 'none')
    drawnow();
    saveas(h_f, [file_name, '_trial.jpeg']);
    saveas(h_f,[file_name, '_trial']);



    %% shift everything according to start_here (first 5 frames when the
    % resonant is not ready)

    start_here = 6;

    quad_data = quad_data(start_here : end);

    field_names = fieldnames(hardware_ind);
    for j = 1 : length(field_names)
        hardware_ind.(field_names{j}) = ...
            hardware_ind.(field_names{j}) - start_here + 1;
        out_of_range = hardware_ind.(field_names{j}) < 1;
        hardware_ind.(field_names{j})(out_of_range) = [];
    end

    remove_us_vec = logical(remove_us_vec(start_here : end));
    save([file_name, '_remove_us_vec_trimmed.mat'],...
        'remove_us_vec');

    quad_bin_number = 10;
    [lap_vec_custom, discrete_quad_data, quad_data_norm] = ...
        create_lap_vec_custom(quad_data, quad_bin_number, hardware_ind,...
        remove_us_vec);
    save([file_name, '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
        'lap_vec_custom', 'discrete_quad_data', 'quad_data_norm',...
        'quad_bin_number');

    %%
    quad_data_smooth = smooth(quad_data, 'moving', 5);
    movement = [0; diff(quad_data_smooth)>1];
    velocity_by_quad_data_smooth = [0; diff(quad_data_smooth)];
    velocity_raw = diff(quad_data);

    quad_data_shift = [];
    save([file_name, '_quad_data_shift_movement_velocity_hardware.mat'],...
        'quad_data_shift', 'movement', 'velocity_raw',...
        'velocity_by_quad_data_smooth', 'hardware_ind',...
        'start_here');

    % legacy : quad_data_shift, velocity_raw;
end



