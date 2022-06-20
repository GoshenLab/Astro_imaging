
mouse_dir = 'D:\astro_imaging\Nature_code\all_code_data\';

template_mat_dir = ...
    'D:\astro_imaging\Nature_code\all_code_data\template_mat_location\';

file_names = {...
    '9N3_000_000',...
    '9O2_000_003',...
    };

only_matching = true;
if only_matching
    name_includes = 'matching_ROIs';
else
    name_includes = 'all_ROIs';
end
ALL_DATA = struct();
ALL_DATA_ind = 1;

for i = 1 : length(file_names)
    mouse_name = file_names{i}(1:3);
    load([mouse_dir, file_names{i} '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'])
    load([mouse_dir, file_names{i} '_quad_data_shift_movement_velocity_hardware.mat'])
    IR_ind = hardware_ind.IR_ind;
    valve_ind = hardware_ind.valve_ind;
    
    [~, weird_lap_number] = max(diff(quad_data_norm(valve_ind)));
    weird_lap_number = weird_lap_number + 1;
    
    lap_type_vec = nan(size(quad_data_norm));
    
    weird_lap_indices = ...
        valve_ind(weird_lap_number) : valve_ind(weird_lap_number + 1);
    familiar_laps_indices = 1 : (weird_lap_indices(1) - 1);
    novel_laps_indices = (weird_lap_indices(end) + 1) : length(quad_data_norm);
    
    lap_type_vec(familiar_laps_indices) = 1;
    lap_type_vec(weird_lap_indices) = 2;
    lap_type_vec(novel_laps_indices) = 3;
    
    figure;
    sgtitle(mouse_name);
    
    plot(quad_data_norm, '.');
    hold on;
    plot(familiar_laps_indices, quad_data_norm(familiar_laps_indices), '.')
    plot(weird_lap_indices, quad_data_norm(weird_lap_indices), '.')
    plot(novel_laps_indices, quad_data_norm(novel_laps_indices), '.')
    plot(IR_ind, quad_data_norm(IR_ind), '^')
    plot(valve_ind, quad_data_norm(valve_ind), 'X')
    legend('path', 'fam', 'weird', 'novel', 'IR', 'valve')
    
    if only_matching 
        template_mat_file = [template_mat_dir, mouse_name, '_' name_includes...
            '_lap_max21_discrete_day_1_vs_day_2_template_mat.mat'];
        template_mat_perm_file = [template_mat_dir, mouse_name, '_' name_includes...
            '_lap_max21_discrete_day_1_vs_day_2_template_mat_perm.mat'];

    else
        template_mat_file = [template_mat_dir, mouse_name, '_' name_includes...
            '_lap_max21_discrete_', file_names{i},...
            '_template_mat.mat'];
        template_mat_perm_file = [template_mat_dir, mouse_name, '_' name_includes...
            '_lap_max21_discrete_', file_names{i},...
            '_template_mat_perm.mat'];
    end
    
    load(template_mat_file)
    load(template_mat_perm_file)
    
    
    PERM_NUM = size(template_mat_perm, 3);
    
    
    h_f = figure();
    sgtitle(mouse_name);

    subplot_ind = 1;
    include_laps = 1 : weird_lap_number - 1;
    disp(file_names{i});
    ALL_DATA = analyse_template_mat_in_parts(...
        template_mat, template_mat_perm,...
        include_laps, subplot_ind,...
        ALL_DATA, ALL_DATA_ind);
    ALL_DATA(ALL_DATA_ind).file_name = file_names{i};
    
    ALL_DATA_ind = ALL_DATA_ind + 1;
    subplot_ind = 1;

    include_laps = ...
        (weird_lap_number + 1) : (weird_lap_number + 10);

    ALL_DATA = analyse_template_mat_in_parts(...
        template_mat, template_mat_perm,...
        include_laps, subplot_ind,...
        ALL_DATA, ALL_DATA_ind);
    ALL_DATA(ALL_DATA_ind).file_name = file_names{i};
    
    ALL_DATA_ind = ALL_DATA_ind + 1;
    
    
    
    
end







function ALL_DATA = analyse_template_mat_in_parts(template_mat, template_mat_perm,...
    lap_ind, plot_ind, ALL_DATA, ALL_DATA_ind)

curr_data_real = template_mat(2:10, lap_ind);
curr_data_perm = template_mat_perm(2:10, lap_ind, :);
PERM_NUM = size(template_mat_perm, 3);

mean_for_each_perm = ...
    squeeze(mean(curr_data_perm, 2, 'omitnan'));
mean_across_all_perm_V1 = mean(mean_for_each_perm, 2, 'omitnan');

subplot(2,1,plot_ind)
shadedErrorBar(2:10, mean(curr_data_real, 2, 'omitnan'),...
    std(curr_data_real,[], 2, 'omitnan')./...
    sqrt(sum(~isnan(curr_data_real), 2)))

subplot(2,1,plot_ind+1)
x = 1 : 9;
shadedErrorBar(x,...
    mean(curr_data_real, 2, 'omitnan') ./ mean_across_all_perm_V1,...
    std(curr_data_real, [], 2, 'omitnan')./...
    sqrt(sum(~isnan(curr_data_real),2))./...
    mean_across_all_perm_V1);
hold on;
plot([min(x) max(x)], [1 1], 'k--')


filler = repmat(2 : 10, size(curr_data_real, 2), 1)';

R_real = corr(filler(:), curr_data_real(:), 'rows', 'pairwise');
R_perm = nan(size(curr_data_perm, 3), 1);
for s = 1 : size(curr_data_perm, 3)
    curr_data = curr_data_perm(:, :, s);
    curr_data = curr_data(:);
    R_perm(s) = corr(filler(:), curr_data(:), 'rows', 'pairwise');
end

ALL_DATA(ALL_DATA_ind).R_real = R_real;
ALL_DATA(ALL_DATA_ind).R_perm = R_perm;

ALL_DATA(ALL_DATA_ind).p_val_new = ...
    sum(R_real < R_perm)./PERM_NUM;
end