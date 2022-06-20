%% run this after you have the 3d_mat!

sim_num = 1000;
lap_max = [];

file_names = {...
    '8C3_100_001';...
    '9B2_100_004';...
    '9H3_003_000';...
    '9I3_004_000';...
    '9P2_000_000';...
    '9N1_000_000';...
    '9O3_000_000';...
    '9Q4_000_000';...
    '9S3_000_000';...
    };

data_folder_name = 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_folder_name);
%%


sig_SI_ROI_num = nan(size(file_names));
sig_SI_ROI_prop = sig_SI_ROI_num;
activation_peak = cell(size(sig_SI_ROI_num));
tuning_curves = activation_peak;


sig_SI_ROI_num_matching = nan(size(file_names));
sig_SI_ROI_prop_matching = sig_SI_ROI_num;
activation_peak_matching = cell(size(sig_SI_ROI_num));
tuning_curves_matching = activation_peak;


for discrete_analysis = [false]
    for j = 1 : length(file_names)
        
        file_name = file_names{j};
        disp(file_name);

        if discrete_analysis
            save_name_3D = [file_name, '_discrete'];
        else
            save_name_3D = [file_name, '_analoge'];
        end
        
        if ~isempty(lap_max)
            save_name = [save_name_3D, '_lap_max_' num2str(lap_max)];
        else
            save_name = save_name_3D;
        end
        
        disp(save_name);
        if ~isfile([save_name, '_sig_SI_ID_final.mat'])
            disp('calculating SI var')
            load([file_name, '_full_sig_ROIs_df_F.signals'], '-mat')
            load([file_name, '_new_event_det_df_F.mat'], 'events_above_min');
            event_mat = double(events_above_min);
            how_many_cells = size(events_above_min, 2);
            load([file_name, '_remove_us_vec_trimmed.mat'])
            
            remove_us_vec = logical(remove_us_vec);
            df_F(remove_us_vec, :) = nan;
            event_mat(remove_us_vec, :) = nan;
            
            if discrete_analysis
                sig = event_mat;
            else
                sig = df_F;
            end
            
            all_ROIs = 1 : how_many_cells;
            
            load([file_name,...
                '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
                'lap_vec_custom',...
                'discrete_quad_data');
            lap_vec = lap_vec_custom;
            load([file_name,...
                '_quad_data_shift_movement_velocity_hardware.mat'],...
                'movement');
            
            movement = logical(movement);
            
            good_frames = ...
                (sum(isnan(sig), 2) == 0) & (movement == 1) &...
                ~isnan(discrete_quad_data');
            
            if ~isempty(lap_max)
                good_frames(lap_vec > lap_max) = 0;
            end
            %%
            if isfile([file_name, '_remove_ROIs.mat'])
                load([file_name, '_remove_ROIs.mat'])
            else
                remove_ROIs = [];
            end
            
            if isfile([file_name, '_elimination_movie_based.mat'])
                load([file_name, '_elimination_movie_based.mat'])
            else
                elimination_movie_based.dead_ROIs = [];
            end
            remove_ROIs = union(remove_ROIs, elimination_movie_based.dead_ROIs);
            
            cell_IDs_3D_mat = setdiff(all_ROIs, remove_ROIs)';
            
            sig = sig(good_frames, cell_IDs_3D_mat);
            
            discrete_quad_data = discrete_quad_data(good_frames);
            
            
            time_spent_in_bin = ...
                splitapply(@length, discrete_quad_data, discrete_quad_data);
            
            time_spent_in_bin = time_spent_in_bin(2:end);
            
            total_time = sum(time_spent_in_bin);
            
            p_bin = [time_spent_in_bin./total_time];
            
            event_prob_in_bin = ...
                splitapply(@mean, sig, discrete_quad_data');
            event_prob_in_bin = event_prob_in_bin(2:end, :);
            
            mean_event_prob = p_bin * event_prob_in_bin;
            
            norm_event_prob = event_prob_in_bin./mean_event_prob;
            norm_event_prob(isequal(norm_event_prob, inf)) = nan;
            SI_real = ...
                sum(p_bin' .* norm_event_prob ...
                .* log2(norm_event_prob), 'omitnan');
            
            SI_nan_real = sum(isnan(SI_real));
            
            perm_delta = randperm(length(sig), sim_num);
            SI_perm = nan(sim_num, size(sig, 2));
            for sim = 1 : sim_num
                if mod(sim, 100) == 0
                    disp(['sim num: '  num2str(sim)])
                end
                curr_sig = circshift(sig, perm_delta(sim));
                event_prob_in_bin = ...
                    splitapply(@mean, curr_sig, discrete_quad_data');
                event_prob_in_bin = event_prob_in_bin(2: end, :);
                mean_event_prob = p_bin * event_prob_in_bin;
                norm_event_prob_perm = event_prob_in_bin./mean_event_prob;
                norm_event_prob_perm(isequal(norm_event_prob_perm, inf)) = nan;
                SI_perm(sim, :) = ...
                    sum(p_bin' .* norm_event_prob_perm ...
                    .* log2(norm_event_prob_perm), 'omitnan');
                
            end
                        
            save([save_name, '_sig_SI_ID_final.mat'],...
                'SI_real', 'SI_perm', 'event_prob_in_bin', 'sim_num',...
                'cell_IDs_3D_mat')
            
        else
            disp('loading SI var')
            load([save_name, '_sig_SI_ID_final.mat'])
            
        end
        
        sig_SI = SI_real > prctile(SI_perm, 95);
        
        load([save_name_3D '_sig_3D_mat_df_F_custom_lap.mat'], 'sig_3D_mat');
        
        curr_cells = cell_IDs_3D_mat(sig_SI);
        if isempty(lap_max)
            curr_3D_data = sig_3D_mat(2:10, :, curr_cells);
        else
            curr_3D_data = ...
                sig_3D_mat(...
                2:10,...
                1 : min(lap_max, size(sig_3D_mat, 2)),...
                curr_cells);
        end
        tuning_curves{j} = squeeze(mean(curr_3D_data, 2, 'omitnan'));
        [~, activation_peak{j}] = max(tuning_curves{j});
    end
end


%% heatmap
figure;
tuning_curves_all = cell2mat(tuning_curves');
tuning_curves_all = sortrows(tuning_curves_all', 5);
subplot(1,2,1)
imagesc(tuning_curves_all)
subplot(1,2,2)
imagesc(normalize(tuning_curves_all,2,'range', [0 1]))
%% bar graph
figure;
A = cell2mat(activation_peak');
h_b = bar(unique(A), histc(A, unique(A)));
xlabel('max event probability location (bin)')
ylabel('ROI number')
h_b.Parent.XTick = [1:9];
h_b.Parent.TickDir = 'out';
h_b.Parent.Box = 'off';

