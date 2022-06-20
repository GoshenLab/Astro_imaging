function mask_struct = create_multi_day_arrays(day_compare_struct, day_fields)

temp = struct2cell(day_compare_struct);

days_with_data = cellfun(@(x) ~isempty(x), temp);
% this analysis is only relevant when you have more than 1 day, so
% there should be at least 3 fields (name + 2 days):
counter = 0; % this is the bias for each day!
mask_struct = struct();
days_with_data = find(days_with_data);
mask_struct.how_many_cells = zeros(length(days_with_data), 1);
sig_array = cell(length(days_with_data), 1);
event_array = sig_array;
movement_array = sig_array; % at the moment, this is based on the "simple" movement calculation

for i = 1 : length(days_with_data)
    
    
    curr_file = day_compare_struct.(day_fields{days_with_data(i)});
    

    %%  not used (12.10.21)
%     load([curr_file '_full_sig_ROIs_df_F.signals'], '-mat');    
% 
%     load([curr_file, '_new_event_det_df_F.mat']);
%     event_mat = events_above_min;


%% this was commented:
%     load([curr_file, '_remove_us_vec_trimmed.mat'])
%     deltaF_dombeck_sig_df_f(remove_us_vec, :) = nan;
%     sig_array{i} = deltaF_dombeck_sig_df_f;

%     load([curr_file, '_quad_data_shift_movement_velocity_hardware.mat'])
%     
%     movement(remove_us_vec) = nan;
%     movement_array{i} = movement;

    %%  not used (12.10.21)

%     df_F(remove_us_vec, :) = nan;
%     sig_array{i} = df_F;
%     event_mat(remove_us_vec, :) = nan;
%     event_array{i} = event_mat;
    

    
    
    
    
    %% now prepare the mask_struct with the required info:
    
    % load the ROI masks for each ot:
    s1 = dir([curr_file '_ot_000_*.segment']);
    s1 = load(s1.name, '-mat', 'mask');
    mask_struct(i).ot_000_ROI_num = max(s1.mask(:));
    
    s1.mask(s1.mask>0) = s1.mask(s1.mask>0) ...
        + counter;
    mask_struct(i).ot_000 = s1.mask;
    
    s2 = dir([curr_file '_ot_001_*.segment']);
    s2 = load(s2.name, '-mat', 'mask');
    mask_struct(i).ot_001_ROI_num = max(s2.mask(:));
    
    % increase the ot number:
    s2.mask(s2.mask>0) = s2.mask(s2.mask>0) ...
        + max(s1.mask(:));
    mask_struct(i).ot_001 = s2.mask;
    
    % add name
    mask_struct(i).name = curr_file;
    
    % count the number of cells in both ots
    mask_struct(i).how_many_cells =  max(s2.mask(:)) -  counter;
    
    % increase the bias, so that the ROI numbers are unique across
    % days:
    counter = max(s2.mask(:));
    
    % 12.4.20: the ROI count in field - should match the full_sig mat
    % (how many cells in each plane...).
    
    %%  load the bad ROIs from the 2 relevant files, and merge into one:
    if isfile([curr_file, '_remove_ROIs.mat'])
    load([curr_file, '_remove_ROIs.mat']);
    else
        remove_ROIs = [];
    end
    
    elimination_file = [curr_file, '_elimination_movie_based.mat'];
    
    % I added this because when I had "perfect" movies, I didn't always
    % save this variable:
    if isfile(elimination_file)
        load(elimination_file);
    else
        elimination_movie_based.dead_ROIs = [];
    end
    remove_ROIs = union(remove_ROIs,...
        elimination_movie_based.dead_ROIs);
    % align to the same orientation
    mask_struct(i).remove_ROIs =...
        reshape(remove_ROIs, length(remove_ROIs), 1);
    
end

save([curr_file(1:3), '_mask_struct_df_F.mat'], 'mask_struct');
%% not saving these because I don't use them (12.10.21)
% save([curr_file(1:3), '_sig_array_df_F.mat'], 'sig_array');
% save([curr_file(1:3), '_event_array_df_F.mat'], 'event_array');
% save([curr_file(1:3), '_movement_array_df_F.mat'], 'movement_array');

end