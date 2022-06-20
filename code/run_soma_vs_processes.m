clear
clc
close all


%% single ROI stuff:

STATs = struct();

file_names = {...
    '8C3_100_001',...
    '9B2_100_004',...
    '9H3_003_000',...
    '9P2_000_000',...
    '9N1_000_000',...
    '9O3_000_000',...
    '9Q4_000_000',...
    '9S3_000_000',...
    };

%    '9I3_004_000',... not included because there aren't enough somas

data_folder  = 'D:\astro_imaging\Nature_code\all_code_data\';
    cd(data_folder)

template_mat_folder = ...
    'D:\astro_imaging\Nature_code\all_code_data\template_mat_location\';

SAMPLING_RATE = 15.49;
PERM_NUM = 1000;
get_event_prob = true;
plot_normalized_ramps = true;
discrete_analysis = true;
bin_num = 10;

for i = 1 : length(file_names)
    
    file_name = file_names{i};
    disp(file_name)
    STATs(i).file_name = file_name;
    
    load([file_name,...
        '_quad_data_shift_movement_velocity_hardware.mat'],...
        'movement');
    load([file_name,...
        '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
        'lap_vec_custom', 'discrete_quad_data', 'quad_data_norm');
    lap_vec = lap_vec_custom;
    movement = logical(movement);
    
    load([file_names{i}, '_new_event_det_df_F.mat'], 'events_above_min');
    event_mat = double(events_above_min);
    load([file_names{i}, '_full_sig_ROIs_df_F.signals'], '-mat',...
        'deltaF_dombeckMat_param');

    load([file_names{i} '_remove_us_vec_trimmed'])
    remove_us_vec = logical(remove_us_vec);
    event_mat(remove_us_vec, :) = nan;
    
    
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
    remove_ROIs = ...
        union(remove_ROIs, elimination_movie_based.dead_ROIs);
    
    
    only_processes_files = dir([file_name, '*only_processes.mat']);
    shift_ot = deltaF_dombeckMat_param.ot_length(1);
    only_processes_log = zeros(sum(deltaF_dombeckMat_param.ot_length), 1);
    for j = 1 : 2
        load(only_processes_files(j).name);
        if j == 2
            only_processes = only_processes + shift_ot;
        end
        only_processes_log(only_processes) = true;
    end
    
    only_processes_cell_IDs = ...
        setdiff(find(only_processes_log), remove_ROIs);
    only_somata_cell_IDs = ...
        setdiff(find(~only_processes_log), remove_ROIs);
    all_cell_IDs_sanity = ...
        setdiff(1 : length(only_processes_log), remove_ROIs);
    
    STATs(i).only_processes = length(only_processes_cell_IDs);
    STATs(i).only_somata = length(only_somata_cell_IDs);
    STATs(i).ROI_num = length(all_cell_IDs_sanity);
    
    
    %%
    discrete_data = discrete_quad_data;
    sig = event_mat;

    %%
    for processes = [true, false] % perform analysis on processes/somata only
        
        if processes == true
            curr_sig = sig(:, only_processes_cell_IDs);
            fig_title = [file_name, '_only_processes'];
        else
            curr_sig = sig(:, only_somata_cell_IDs);
            fig_title = [file_name, '_only_somata'];
        end
        
        concurrent_event_num = sum(curr_sig, 2);        

        good_frame_ind = ...
            ~((sum(isnan(curr_sig), 2) > 0) |...
            (movement ~= 1));
        
        template_mat = ...
            number_of_sim_events_in_bin...
            (lap_vec(good_frame_ind),...
            discrete_data(good_frame_ind),...
            concurrent_event_num(good_frame_ind), 1);
        
        template_mat = template_mat(2:10, :, :);

        % this is in order to get the cell proportion!
        template_mat = template_mat./size(curr_sig, 2);

        
        if ~isfile([template_mat_folder, fig_title, '_template_mat_perm.mat'])
            disp('creating template mat perm')
            sig_for_perm = curr_sig(good_frame_ind, :);
            perm_delta = randperm(length(sig_for_perm), PERM_NUM);
            template_mat_perm = ...
                nan(bin_num, max(lap_vec), PERM_NUM);
            for p = 1 : PERM_NUM
                
                concurrent_event_num_perm = ...
                    circshift(concurrent_event_num(good_frame_ind),...
                    perm_delta(p));
                
                temp = number_of_sim_events_in_bin...
                    (lap_vec(good_frame_ind),...
                    discrete_data(good_frame_ind),...
                    concurrent_event_num_perm, 1);
                
                %% normalization for discrete analysis:
                temp = temp./size(curr_sig, 2);
                %%
                template_mat_perm(1 : size(temp, 1),...
                    1 : size(temp, 2), p) = ...
                    temp;
            end
            save([template_mat_folder, fig_title, '_template_mat_perm.mat'],...
                'template_mat_perm', 'PERM_NUM')
        else
            load([template_mat_folder, fig_title, '_template_mat_perm.mat'],...
                'template_mat_perm', 'PERM_NUM')
            disp('loading template mat perm')
            
        end
            
        template_mat_perm = template_mat_perm(2:10, :, :);

        mean_for_each_perm = ...
            squeeze(mean(template_mat_perm, 2, 'omitnan'));
        mean_across_all_perm_V1 = mean(mean_for_each_perm, 2, 'omitnan');

        if plot_normalized_ramps
           
            x = [2 : 10];
            figure;
            shadedErrorBar(x,...
                mean(template_mat, 2, 'omitnan') ./ mean_across_all_perm_V1,...
                std(template_mat, [], 2, 'omitnan')./...
                sqrt(sum(~isnan(template_mat),2))./...
                mean_across_all_perm_V1);
            hold on;
            plot([min(x) max(x)], [1 1], 'k--')
            title({[fig_title, ' normalized by permutation average']},...
                'interpreter', 'none');
            
        end
        
        curr_data_real = template_mat;
        curr_data_perm = template_mat_perm;
        filler = repmat(2 : 10, size(template_mat, 2), 1)';
    
        R_real = corr(filler(:), curr_data_real(:), 'rows', 'pairwise');
        filler = repmat(2 : 10, size(curr_data_perm, 2), 1)';
        R_perm = nan(size(curr_data_perm, 3), 1);
        for s = 1 : size(curr_data_perm, 3)
            curr_data = curr_data_perm(:, :, s);
            curr_data = curr_data(:);
            R_perm(s) = corr(filler(:), curr_data(:), 'rows', 'pairwise');
        end
    
        if processes
            struct_title = 'only_processes';
        else
            struct_title = 'only_somata';
        end
        STATs(i).([struct_title, '_R_real']) = R_real;
        STATs(i).([struct_title, '_R_real_p_val']) = sum(R_real<R_perm)/PERM_NUM;
    end
    
end

%%
pooled_processes = [STATs.only_processes_R_real];
pooled_somas = [STATs.only_somata_R_real];


h_f = figure;
h_p = plot([ones(size(pooled_somas)); 2*ones(size(pooled_processes))],...
    [pooled_somas; pooled_processes], 'color', [.5 .5 .5]);
hold on;
plot([ones(size(pooled_somas)); 2*ones(size(pooled_processes))],...
    [pooled_somas; pooled_processes], '.', 'markersize', 25)
xlim([0 3])
ylim([0 0.6])
h_f.Children.Box = 'off';
h_f.Children.TickDir = 'out';
