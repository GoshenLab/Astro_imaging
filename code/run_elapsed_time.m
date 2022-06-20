clear
clc
close all


file_names = {...
    '8C3_100_001',...
    '9B2_100_004',...
    '9H3_003_000',...
    '9I3_004_000',...
    '9P2_000_000',...
    '9N1_000_000',...
    '9O3_000_000',...
    '9Q4_000_000',...
    '9S3_000_000',...
    };

data_folder = 'D:\astro_imaging\Nature_code\all_code_data\';


SAMPLING_RATE = 15.49;
PERM_NUM = 1000;

SIG_struct = struct();

samples_after_rew = round(SAMPLING_RATE * 7);
discrete_analysis = [true];

c = 1;

for m = 1 : length(file_names)
    
    file_name = file_names{m};
    save_name = [file_name, '_discrete'];
   
    
    %% first deal with the signal, which remains constant:
    load([file_name, '_new_event_det_df_F.mat'],...
        'events_above_min');
    event_mat = double(events_above_min);
    
    cell_IDs = 1 : size(event_mat, 2);
    
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

    cell_IDs = setdiff(cell_IDs, remove_ROIs);
    
    event_mat = event_mat(:, cell_IDs);
    load([file_name, '_remove_us_vec_trimmed.mat'])
    
    remove_us_vec = logical(remove_us_vec);
    event_mat(remove_us_vec, :) = nan;
    
    
    %% 9.8.20: I use the lap_vec_custom variables instead...
    load([file_name,...
        '_quad_data_shift_movement_velocity_hardware.mat'],...
        'movement', 'hardware_ind');

    movement = logical(movement);
    
    sig = event_mat;

    concurrent_event_num = sum(sig, 2);
    
    reward_ind = ...
        setdiff(hardware_ind.valve_ind, find(remove_us_vec));
    
    %%
    after_rew = nan(length(reward_ind), samples_after_rew);
    for r = 1 : length(reward_ind)
        if (reward_ind(r) + samples_after_rew - 1) <= ...
                length(concurrent_event_num)
            after_rew(r, :) = ...
                concurrent_event_num(reward_ind(r) : ...
                reward_ind(r) + samples_after_rew - 1);
        end
    end
 
    after_rew_prop = after_rew./(length(cell_IDs));
    
    after_rew_prop_perm = nan(size(after_rew_prop, 1),...
        size(after_rew_prop, 2),...
        PERM_NUM);
    
%     if numel(after_rew_prop)< PERM_NUM
%         continue
%     end
    cyc_shift = randperm(numel(after_rew_prop), PERM_NUM);
    
    perm_me = after_rew_prop';
    perm_me = perm_me(:);
    for s = 1 : PERM_NUM
        tempi = circshift(perm_me, cyc_shift(s));
        tempi = ...
            reshape(tempi, size(after_rew_prop, 2), size(after_rew_prop, 1));
        after_rew_prop_perm(:, :, s) = tempi';
    end
    
    filler = repmat(1 : size(after_rew_prop, 2),...
        size(after_rew_prop, 1), 1);
    
    
    R_real = corr(filler(:), after_rew_prop(:), 'rows', 'pairwise');
    R_perm = nan(size(after_rew_prop_perm, 3), 1);
    for s = 1 : size(after_rew_prop_perm, 3)
        curr_data = after_rew_prop_perm(:, :, s);
        curr_data = curr_data(:);
        R_perm(s) = corr(filler(:), curr_data(:), 'rows', 'pairwise');
    end
    
    SIG_struct(c).file_name = file_name;
    SIG_struct(c).R_real = R_real;
    
    SIG_struct(c).p_val_new = sum(R_perm<R_real)./PERM_NUM;
      
    curr_data = after_rew_prop./mean(after_rew_prop(:), 'omitnan');   

    subplot(3, 3, c)
    shadedErrorBar([], mean(curr_data, 'omitnan'),...
        std(curr_data, 'omitnan')./sqrt(sum(~isnan(curr_data))))
    hold on;
    line([0  size(curr_data, 2)-1], [1 1])
    
     c = c + 1;

end



