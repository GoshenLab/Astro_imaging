
dist_bin_um = 20;


file_names = {...
    '8C3_100_001',...
    '9B2_100_004',...
    '9H3_003_000',...
    '9I3_004_000',...
    '9P2_000_000',...
    '9N1_000_000',...
    '9O3_000_000',...
    '9Q4_000_000',...
    };


data_folder = 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_folder);

ROI_dist_pooled = [];
ROI_size_pooled = ROI_dist_pooled;
event_corr_pooled = ROI_dist_pooled;
event_prob_pooled = ROI_dist_pooled;

for i = 1 : length(file_names)
    
    disp(file_names{i})
    load([file_names{i} '.mat'])
    mouse_name = file_names{i}(1:3);
    
    load([mouse_name '_mask_struct_df_F_all_ROIs.mat']);
    curr_ind = find(strcmp({mask_struct.name}, file_names{i}));
    
    ROIs_in_prev = 0;
    if curr_ind > 1
        ROIs_in_prev = sum([mask_struct(1 : curr_ind - 1).how_many_cells]);
    end
    
    shift_ROIs = ...
        mask_struct(curr_ind).ot_000_ROI_num;
    
    load([file_names{i}, '_full_sig_ROIs_df_F.signals'], '-mat')
    load([file_names{i}, '_new_event_det_df_F.mat'], 'events_above_min');
    event_mat = double(events_above_min);
    
    load([file_names{i} '_remove_us_vec_trimmed'])
    remove_us_vec = logical(remove_us_vec);
    event_mat(remove_us_vec, :) = nan;
    
    bad_ROIs = mask_struct(1).remove_ROIs;
    event_mat(:, bad_ROIs) = nan;

    if isfile([file_names{i} '_event_mat_corr.mat'])
        disp('loading events correlations')
        
        load([file_names{i} '_event_mat_corr.mat'])
    else
        
        disp('calculating events correlations')
        event_mat_corr = corr(event_mat, 'Rows', 'pairwise');
        save([file_names{i} '_event_mat_corr.mat'], 'event_mat_corr')
    end
    
    figure;
    sgtitle(file_names{i}, 'interpreter', 'none')
    
    scaled_centroids = nan(size(event_mat, 2), 2);
    ROI_size = nan(size(event_mat, 2), 1);
    
    for ot = 0:1
        curr_mask = mask_struct(curr_ind).(['ot_00', num2str(ot)]);
        
        if ot == 0
            ROI_range = 1 : shift_ROIs;
            curr_mask_copy = zeros(size(curr_mask));
            curr_mask_copy(curr_mask>0) = ...
                curr_mask(curr_mask>0) - (ROIs_in_prev);
            curr_mask = curr_mask_copy;
        else
            ROI_range = (shift_ROIs + 1) : size(df_F, 2);
            curr_mask_copy = zeros(size(curr_mask));
            curr_mask_copy(curr_mask>0) = ...
                curr_mask(curr_mask>0) - (shift_ROIs + ROIs_in_prev);
            curr_mask = curr_mask_copy;
            
        end
        
        
        [~, scaled_centroids(ROI_range, :), FOV_size] = ...
            find_ROI_centroid(curr_mask, info);
        ROI_size_temp = histc(curr_mask(:), unique(curr_mask(:)));
        ROI_size_temp = ROI_size_temp(2:end); % omit the zero BG
        ROI_size(ROI_range) = ROI_size_temp;
    end
    
    event_prob = mean(event_mat, 'omitnan');
    
    subplot(1,2,1)
    scatter(ROI_size, event_prob);
    xlabel('ROI size (pixels)')
    ylabel('event probability')
    
       
    %%
    
    scaled_centroids_dist = pdist2(scaled_centroids, scaled_centroids);
    tri_extract = tril(true(length(scaled_centroids_dist)),-1);
    ot_1_ind = [(mask_struct(curr_ind).ot_000_ROI_num+1), length(event_mat_corr)];
    ot_0_ind = [1 ,(mask_struct(curr_ind).ot_000_ROI_num)];
    tri_extract(ot_1_ind(1): ot_1_ind(2), ot_0_ind(1):ot_0_ind(2)) = 0;

    scaled_centroids_dist = scaled_centroids_dist(tri_extract);
    curr_event_mat_corr = event_mat_corr(tri_extract);

    


ROI_dist_pooled = [ROI_dist_pooled; scaled_centroids_dist];
ROI_size_pooled = [ROI_size_pooled; ROI_size];
event_corr_pooled = [event_corr_pooled; curr_event_mat_corr];
event_prob_pooled = [event_prob_pooled; event_prob'];
end

figure;
h_ax = subplot(1,2,1);
create_graphs_discrete_X(ROI_dist_pooled, 0:10:1000, event_corr_pooled,...
    h_ax, true, 300)
xlabel('ROI centroid distance (um)')
ylabel('mean event correlation')
ylim([0 0.3])
h_ax.Box = 'off';
h_ax.TickDir = 'out';

h_ax = subplot(1,2,2);
create_graphs_discrete_X(ROI_size_pooled, 0:10:2000, event_prob_pooled,...
    h_ax, true, 10)
xlabel(h_ax, 'ROI size (pixels)')
ylabel(h_ax, 'mean event probability')
ylim([0 0.5])
h_ax.Box = 'off';
h_ax.TickDir = 'out';


function create_graphs_discrete_X(dist_vec, bin_vec, ana_var, h_ax,...
    pooled_data, thresh)


dist_binned = discretize(dist_vec, bin_vec);

bin_vec_plot = bin_vec(1 : end-1) + 0.5 * mean(diff(bin_vec));

mean_corr_per_dist = ...
    splitapply(@(x) mean(x, 'omitnan'), ana_var,...
    findgroups(dist_binned));
std_corr_per_dist = ...
    splitapply(@(x) std(x, 'omitnan'), ana_var,...
    findgroups(dist_binned));
N_corr_per_dist = ...
    splitapply(@(x) sum(~isnan(x)), ana_var,...
    findgroups(dist_binned));
bin_num = unique(dist_binned(~isnan(dist_binned)));

if ~pooled_data
    errorbar(h_ax, bin_vec_plot(bin_num),...
        mean_corr_per_dist, std_corr_per_dist./sqrt(N_corr_per_dist),...
        'Color', 'k', 'LineStyle', 'none')
else
    ind = N_corr_per_dist > thresh;
    errorbar(h_ax, bin_vec_plot(bin_num(ind)),...
        mean_corr_per_dist(ind), std_corr_per_dist(ind)./sqrt(N_corr_per_dist(ind)),...
        'Color', 'k', 'LineStyle', 'none')   
end
end
