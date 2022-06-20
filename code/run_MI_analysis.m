close all
clear
clc

%%
data_folder = 'D:\astro_imaging\Nature_code\all_code_data\';

graph_matching_results_folder = ...
    'D:\astro_imaging\Nature_code\all_code_data\matching_days_data\';

mouse_names = {
    '9Q4',...
    '9P2',...
    '8C3',...
    '9B2',...
    '9O3',...
    '9N1',...
    '9S3',...
    };


sim_num = 1000; 
lap_max = 20;

ALL_DATA = struct();
sig_pairs_prop = nan(size(mouse_names, 1), 2);

c = 1; 

for m = 1 : length(mouse_names)
    mouse_name = mouse_names{m};
    
    load([mouse_name '_mask_struct_df_F_matching_ROIs.mat'])
    load([graph_matching_results_folder mouse_name, '_graph_match_mat_df_F.mat']);
    
    cells_per_day = [mask_struct.how_many_cells];
    day_shift = [0, cumsum(cells_per_day)]';
    
    file_names = {mask_struct.name};
    file_names = file_names(1 : 2);
    
    g_active = rmnode(g, find(~ismember(g.Nodes.day_ind, [1 2])));
    g_active = rmnode(g_active, find(g_active.Nodes.remove_me == 1));
    g_active = rmnode(g_active, find(g_active.degree == 0));
    
    matching_ROIs = cellfun(@str2double, g_active.Edges.EndNodes);
    matching_ROIs(:, 1) = matching_ROIs(:, 1) - day_shift(1);
    matching_ROIs(:, 2) = matching_ROIs(:, 2) - day_shift(2);
    
    if strcmp(mouse_name, '9Q4')
        % matching ot_000 only!
        omit_ot = [mask_struct.ot_000_ROI_num];
        A = (matching_ROIs <= omit_ot(1:2));
        matching_ROIs(sum(A, 2) == 0, :) = [];
    end
    
    
    for i = 1 : length(file_names)
        file_name = file_names{i};
        disp(file_name)
        save_name = [file_name, '_MI_analysis_lap_max_' num2str(lap_max)];

        if isfile([save_name '.mat'])
            disp(['loading ' save_name])
            load([save_name])
        else
            disp(['creating ' save_name])
            
            load([file_name '_new_event_det_df_F.mat'])
            load([file_name, '_remove_us_vec_trimmed.mat'])
            load([file_name, '_remove_ROIs.mat'])
            load([file_name, '_quad_data_shift_movement_velocity_hardware.mat'],...
                'movement')
            load([file_name, '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
                'lap_vec_custom')
            
            sig_trim = events_above_min;
            remove_vec = zeros(size(remove_us_vec));
            
            remove_vec(remove_us_vec) = true;
            remove_vec(movement'~=1) = true;

            if ~isempty(lap_max)
               remove_vec(lap_vec_custom > lap_max) = true;
            end
            
            sig_trim(logical(remove_vec), :) = [];
            sig_trim(sum(isnan(sig_trim), 2) > 0, :) = [];
            
            % the match mat is after node elimination!
            % the sig_trim is organized according to ROI id so that the MI
            % matrices are for matching ROIs across days.
            
            sig_trim = sig_trim(:, matching_ROIs(:, i));
            
            MI_mat = calculate_MI_mat(sig_trim);
            MI_mat_perm = create_MI_mat_perm(sig_trim, sim_num);
            
            save([save_name '.mat'],...
                'MI_mat', 'MI_mat_perm');
        end
        MI_perm_95 = prctile(MI_mat_perm, 95, 3);
        sig_large_MI = MI_mat > MI_perm_95;
        
        ALL_DATA(c).name = file_name;
        ALL_DATA(c).possible_pairs = ...
            (size(MI_mat, 1)^2 - size(MI_mat, 1))/2;
        ALL_DATA(c).sig_large_MI_pairs = sum(sig_large_MI, 'all')/2;
        ALL_DATA(c).sig_large_MI_pairs_prop = ...
            ALL_DATA(c).sig_large_MI_pairs/...
            ALL_DATA(c).possible_pairs;
        
        sig_pairs_prop(m, i) = ALL_DATA(c).sig_large_MI_pairs_prop;
        
        c = c + 1;
    end
    

end


%% ploting:
figure;
plot(sig_pairs_prop', '.', 'MarkerSize', 30)
hold on
plot(sig_pairs_prop', 'Color', [.5 .5 .5])
xlim([0 3])
ylim([0 1])
ylabel('proportion of sig pairs')
xticks([1 2])
xticklabels({'familiar', 'novel'})

