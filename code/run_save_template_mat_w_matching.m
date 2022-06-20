
clear
clc
close all;

ALL_DATA = struct();
only_matching = true;
lap_cutoff = 20;    % could be [] if you want all laps.
quad_or_vel = 1;    % 1: location, 2: velocity
random_reward = false;

if only_matching
    graph_matching_results_folder = ...
        'D:\astro_imaging\Nature_code\all_code_data\matching_days_data\';

    %% familiar-novel:
mouse_names = {'9N1'};
    
    %     mouse_names = {
%         '8C3',...
%         '9B2',...
%         '9P2',...
%         '9N1',...
%         '9O3',...
%         '9Q4',...
%         '9S3',...
%         };


%% reward_shift (max_lap = 21)
% % the first session includes before/after change, so the graph here is
% % missleading - use the run_analyse_template_mat_changing_rew script for
% % analysis.
%     mouse_names = {
%         '9O2',...
%         '9N3',...
%         };


else
    %% for familiar-nov comparison
    mouse_names = {...
        '8C3',...
        '9B2',...
        '9H3',...
        '9I3',...
        '9P2',...
        '9N1',...
        '9O3',...
        '9Q4',...
        '9S3',...
        };
    %% for no VR:
%     mouse_names = {...
%         '9P2',...
%         '9Q4',...
%         '9N3',...
%         '9O2',...
%         };

end

if random_reward
    mouse_names = {...
        '9N1',...
        '9N3',...
        '9O2'...
        };
    only_matching = false;
end
mouse_dir = {...
        'D:\astro_imaging\Nature_code\all_code_data\',...
        };

% this is where you save the results:
folder_name = ...
    'D:\astro_imaging\Nature_code\all_code_data\template_mat_location\';

PERM_NUM = 1000;

bin_number_vel = 10;
bin_number_quad = 10;

c = 1;

for m = 1 : length(mouse_names)
    
    mouse_name = mouse_names{m};
    disp(mouse_name)

    if only_matching
        load([mouse_name, '_mask_struct_df_F_matching_ROIs.mat']);
    else
        load([mouse_name, '_mask_struct_df_F_all_ROIs.mat']);
    end
    matching_days = length(mask_struct);
    cells_per_day = [mask_struct.how_many_cells];
    file_names = {mask_struct.name};
    
    for discrete_analysis = [true]
        
        
        if only_matching
            load([graph_matching_results_folder, mouse_name,...
                '_graph_match_mat_df_F.mat'])
            save_name = [folder_name, mouse_name, '_matching_ROIs'];
            
            day_comp = nan(matching_days - 1, 2);
            day_shift = [0, cumsum(cells_per_day)]';
            
            for d = 1 : matching_days - 1
                day_comp(d, :) = [d, d+1];
            end

            if matching_days == 3
                file_names = ...
                    {file_names{1}, file_names{2},...
                    file_names{2}, file_names{3}};
                matching_days = 4;
            end

            match_mat = cell(matching_days, 1);
            helper = match_mat;

            for k = 1 : size(day_comp, 1)
                g_active = rmnode(g, find(~ismember(g.Nodes.day_ind, day_comp(k, :))));
                g_active = rmnode(g_active, find(g_active.Nodes.remove_me == 1));
                g_active = rmnode(g_active, find(g_active.degree == 0));
                
                
                temp = cellfun(@str2double, g_active.Edges.EndNodes);
                temp(:, 1) = temp(:, 1) - day_shift(day_comp(k, 1));
                temp(:, 2) = temp(:, 2) - day_shift(day_comp(k, 2));
                
                if strcmp(mouse_name, '9Q4') && isequal(day_comp(k, :), [1 2])
                    % matching ot_000 only!
                    ROIs_per_day = [mask_struct.ot_000_ROI_num];
                    A = (temp <= ROIs_per_day(day_comp(k, :)));
                    temp(sum(A, 2) == 0, :) = [];
                end
                
                match_mat{k*2-1} = temp(:, 1);
                match_mat{k*2} = temp(:, 2);
                
                
                helper{k*2-1} = ...
                    ['day_' num2str(day_comp(k,1)) ...
                    '_vs_day_' num2str(day_comp(k,2))];
                helper{k*2} = ...
                    ['day_' num2str(day_comp(k,2)) ...
                    '_vs_day_' num2str(day_comp(k,1))];
            end
            
            
        else
            match_mat = cell(matching_days, 1);
            for i = 1 : matching_days
                match_mat{i} = 1 : cells_per_day(i);
            end
            save_name = [folder_name, mouse_name, '_all_ROIs'];
        end
        
        if ~isnan(lap_cutoff)
            save_name = ...
                [save_name '_lap_max' num2str(lap_cutoff)];
        end
        
        disp(save_name)
        
        if discrete_analysis
            save_name = [save_name, '_discrete'];
        else
            save_name = [save_name, '_analoge'];
        end
        
        %%
        
        for j = 1 : matching_days % change this according to the relevant sessions
           file_name = file_names{j};
            disp(file_name);         
            if random_reward
                if sum(strcmp(file_name, ...
                        {'9N1_003_000', '9N3_003_000', '9O2_003_003'}))==0
                    continue
                end
            end
            if only_matching == 0
                save_name_curr = [save_name, '_', file_name];
            else
                save_name_curr = [save_name, '_', helper{j}];
            end
            if quad_or_vel == 2
                save_name_curr = [save_name_curr '_velocity'];
            end
            
            include_cells = match_mat{j};
            if ~isfile([save_name_curr, '_template_mat.mat']) ||...
                    ~isfile([save_name_curr, '_template_mat_perm.mat'])
                
                [template_mat, template_mat_perm] = ...
                    deal_with_file(file_name, include_cells, lap_cutoff,...
                    quad_or_vel, discrete_analysis, PERM_NUM,...
                    bin_number_vel, bin_number_quad);
                % template mat in the discrete case: prop!!
                
                disp('saving template_mat')
                save([save_name_curr, '_template_mat.mat'],...
                    'template_mat')
                
                disp('saving template_mat_perm')
                save([save_name_curr, '_template_mat_perm.mat'],...
                    'template_mat_perm', 'PERM_NUM')
            else
                disp('loading template_mat')
                load([save_name_curr, '_template_mat.mat'],...
                    'template_mat')
                
                disp('loading template_mat_perm')
                load([save_name_curr, '_template_mat_perm.mat']);
            end
            
            
            % to get the number of cells, not the proportion (not really relevant, unless you want to plot it):
            if discrete_analysis  % include_cells doesn't take into consideration remove_ROIs
                if only_matching == 0
                    template_mat = ...
                        template_mat.* (cells_per_day(j) - ...
                        length(mask_struct(j).remove_ROIs));
                    
                    template_mat_perm = ...
                        template_mat_perm.* (cells_per_day(j) - ...
                        length(mask_struct(j).remove_ROIs));
                else % include_cells DOES take into consideration remove_ROIs
                    template_mat = ...
                        template_mat.* length(include_cells);
                    
                    template_mat_perm = ...
                        template_mat_perm.* length(include_cells);
                    
                end
            end
            
            mean_for_each_perm = ...
                squeeze(mean(template_mat_perm, 2, 'omitnan'));
            mean_across_all_perm_V1 = mean(mean_for_each_perm, 2, 'omitnan');
            

            
            %%
            if  quad_or_vel == 1 && ~random_reward
                x = 2 : 10;
            elseif  quad_or_vel == 1 && random_reward
                x = [2 : 10, 1];    % such that the reward is in bin=10
            elseif quad_or_vel == 2
                x = 1 : 10;
            end
            template_mat = template_mat(x, :);
            mean_across_all_perm_V1 = mean_across_all_perm_V1(x);

            for_nature = template_mat./mean_across_all_perm_V1;

            subplot(3,3,c)
            [~,name_for_title,~] = fileparts(save_name_curr);
            plot_template_mat_results(template_mat, mean_across_all_perm_V1,...
                {file_name; name_for_title})

            filler = repmat(1 : length(x), size(template_mat, 2), 1)';
            
            R_real = corr(filler(:), template_mat(:), 'rows', 'pairwise');
            R_perm = nan(size(template_mat_perm, 3), 1);
            for s = 1 : size(template_mat_perm, 3)
                curr_data = template_mat_perm(x, :, s);
                curr_data = curr_data(:);
                R_perm(s) = corr(filler(:), curr_data(:), 'rows', 'pairwise');
            end
            
            ALL_DATA(c).file_name = file_names{j};
            ALL_DATA(c).only_matching = only_matching;
            ALL_DATA(c).R_real = R_real;
            ALL_DATA(c).r_perm_95_prctile = prctile(R_perm, 95);

            
            ALL_DATA(c).r_real_est_new_p_val = ...
                sum(R_perm > R_real)./PERM_NUM;
            
            c = c + 1;
            
        end
    end
    
end






function [template_mat, template_mat_perm] = ...
    deal_with_file(file_name, include_ROIs, lap_cutoff, quad_or_vel,...
    discrete_analysis, PERM_NUM, bin_number_vel, bin_number_quad)


%% first deal with the signal, which remains constant:
load([file_name, '_full_sig_ROIs_df_F.signals'], '-mat')
load([file_name, '_new_event_det_df_F.mat'],...
    'events_above_min');
event_mat = double(events_above_min);

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

cell_IDs = setdiff(include_ROIs, remove_ROIs); % redundant for the matching case
df_F = df_F(:, cell_IDs);

event_mat = event_mat(:, cell_IDs);
load([file_name, '_remove_us_vec_trimmed.mat'])

remove_us_vec = logical(remove_us_vec);
df_F(remove_us_vec, :) = nan;
event_mat(remove_us_vec, :) = nan;


load([file_name,...
    '_quad_data_shift_movement_velocity_hardware.mat'],...
    'movement', 'velocity_by_quad_data_smooth');
load([file_name,...
    '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
    'lap_vec_custom', 'discrete_quad_data');
lap_vec = lap_vec_custom;

if ~isnan(lap_cutoff)
    cut_point = find(lap_vec == (lap_cutoff + 1), 1);
    if ~isempty(cut_point)
        df_F(cut_point:end, :) = nan;
        event_mat(cut_point:end, :) = nan;
        % to avoid non-existing values
        velocity_by_quad_data_smooth(cut_point:end, :) = nan;
    end
end

movement = logical(movement);
velocity_by_quad_data_smooth(~movement) = nan;


velocity_by_quad_data_smooth(...
    isnan(discrete_quad_data)) = nan; % to include the remove_us_vec
vel_data_norm = ...
    normalize(velocity_by_quad_data_smooth, 'range');
discrete_vel_data = ...
    discretize(vel_data_norm, bin_number_vel);

if quad_or_vel == 1
    discrete_data = discrete_quad_data;

    bin_num = bin_number_quad;

else
    discrete_data = discrete_vel_data;
    bin_num = bin_number_vel;
end

if discrete_analysis
    sig = event_mat;
else
    sig = df_F;
end
%%

if discrete_analysis
    concurrent_event_num = sum(sig, 2);
else
    concurrent_event_num = mean(sig, 2);
end

good_frame_ind = ...
    ~((sum(isnan(sig), 2) > 0) |...
    (movement ~= 1));


template_mat = ...
    number_of_sim_events_in_bin...
    (lap_vec(good_frame_ind),...
    discrete_data(good_frame_ind),...
    concurrent_event_num(good_frame_ind), 1);


if discrete_analysis
    % this is in order to get the cell proportion!
    template_mat = template_mat./size(sig, 2);
    
end

%% permutation
sig_for_perm = sig(good_frame_ind, :);

perm_delta = randperm(length(sig_for_perm), min(PERM_NUM, size(sig_for_perm, 1)));
disp('calculating template_mat_perm')
template_mat_perm = ...
    nan(bin_num, size(template_mat, 2), min(PERM_NUM, size(sig_for_perm, 1)));
for p = 1 : min(PERM_NUM, size(sig_for_perm, 1))
    
    concurrent_event_num_perm = ...
        circshift(concurrent_event_num(good_frame_ind),...
        perm_delta(p));
    
    temp = number_of_sim_events_in_bin...
        (lap_vec(good_frame_ind),...
        discrete_data(good_frame_ind),...
        concurrent_event_num_perm, 1);
    %
    %% normalization for discrete analysis:
    if discrete_analysis
        temp = temp./size(sig, 2);
    end
    %
    template_mat_perm(1 : size(temp, 1),...
        1 : size(temp, 2), p) = ...
        temp;
end


end



function plot_template_mat_results(template_mat, mean_across_all_perm_V1,...
    file_name)

% figure;
x = [1 : size(template_mat, 1)];
shadedErrorBar(x,...
    mean(template_mat, 2, 'omitnan') ./ mean_across_all_perm_V1,...
    std(template_mat, [], 2, 'omitnan')./...
    sqrt(sum(~isnan(template_mat),2))./...
    mean_across_all_perm_V1);
hold on;
plot([min(x) max(x)], [1 1], 'k--')
file_name{end+1} = 'normalized by permutation average';
title(file_name,...
    'interpreter', 'none');
end