clear
clc
% close all
%%

data_folder = 'D:\astro_imaging\Nature_code\all_code_data\';
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


GLM_all = struct();
GLM_all_shuffle = struct();

do_cross_val = false;
GLM_cross_val_results = struct();
GLM_cross_val_shuffle_results = struct();

partial_model_predictor_ind = [2:5]; % column1 = location
perm_num = 1000;

for f = 1 : length(file_names)
    file_name = file_names{f};

    GLM_all(f).name = file_name;
    GLM_cross_val_results(f).name = file_name;

    load([file_name, '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'])
    load([file_name, '_quad_data_shift_movement_velocity_hardware.mat'],...
        'velocity_by_quad_data_smooth')
    load([file_name '_new_event_det_df_F.mat'], 'events_above_min')
    load([file_name '_quad_data_shift_movement_velocity_hardware.mat'],...
        'hardware_ind')
    load([file_name, '_remove_ROIs.mat'])
    if isfile([file_name, '_elimination_movie_based.mat'])
        load([file_name, '_elimination_movie_based.mat'])
        remove_ROIs = union(elimination_movie_based.dead_ROIs, remove_ROIs);
    end
    load([file_name, '_quad_data_shift_movement_velocity_hardware.mat'],...
        'movement')

    %%
    total_lap_number = max(lap_vec_custom);
    elapsed_time = nan(size(lap_vec_custom)); % for now: in samples (starting at zero)
    for i = 1 : total_lap_number
        curr_ind = find(lap_vec_custom == i);
        elapsed_time(curr_ind : end) = ...
            0 : length(curr_ind : length(elapsed_time) - 1);
    end

    %%
    win_length = 30;
    integration_filter = ones(win_length, 1);
    velocity_integral = ...
        conv2(velocity_by_quad_data_smooth, integration_filter, 'same');

    %%
    licking = hardware_ind.lick_ind;
    licking_vec = zeros(size(lap_vec_custom));
    licking_vec(licking) = 1;

    %%
    sig = events_above_min;
    sig(:, remove_ROIs) = [];
    concurrent_events = sum(sig, 2);

    good_sig = ~isnan(concurrent_events); % for bad registration removal

    %%

    X = [...
        quad_data_norm',...
        elapsed_time',...
        velocity_by_quad_data_smooth,...
        velocity_integral,...
        licking_vec'...
        ];

    include_frames = ...
        find((sum(isnan(X),2) == 0) & (movement == 1) & good_sig);
%     if strcmp('8C3_100_001', file_name)
%         include_frames = find((sum(isnan(X),2) == 0) & (movement == 1) & ...
%             good_sig & ~(lap_vec_custom' == 62));
%     end

    full_dataset = X(include_frames, :);
    Y = concurrent_events(include_frames);
    N_samples = length(full_dataset);

    %% testing on all the dataset, wo cross val:

    disp(file_name);
    [GLM_all(f).R_square_full_model,...
        GLM_all(f).R_square_wo_location_model,...
        GLM_all(f).coefficient_p_val,...
        GLM_all(f).RMSE_full_model,...
        GLM_all(f).RMSE_wo_location_model,...
        GLM_all(f).LLR] = ...
        linear_prediction(full_dataset, Y, partial_model_predictor_ind);

    % shuffle the location vector:
    perm_ind = randi(N_samples, [perm_num, 1]);

    for p = 1 : perm_num
        if mod(p, 100) == 0
            disp(p)
        end
        curr_full_dataset = full_dataset;
        curr_full_dataset(:, 1) = ...            % LOCATION is in the first column!
            circshift(full_dataset(:, 1), perm_ind(p));
        [GLM_all_shuffle.(['file_' file_name])(p).R_square_full_model,...
            ~,...                                % I don't permute any input, so it's redundant
            GLM_all_shuffle.(['file_' file_name])(p).coefficient_p_val,...
            GLM_all_shuffle.(['file_' file_name])(p).RMSE_full_model,...
            GLM_all_shuffle.(['file_' file_name])(p).RMSE_wo_location_model] = ...
            linear_prediction(curr_full_dataset, Y, partial_model_predictor_ind);

    end


    %% cross-validation
    disp(file_name);
    N_blocks = 10;

    if do_cross_val
        block_size = 31; % 2 sec
        log_vec_perm = ...
            create_log_vec_uniform_with_random_assignment(...
            N_blocks, N_samples, block_size, perm_num);


        for p = 1 : perm_num
            if mod(p, 100) == 0
                disp(p)
            end

            curr_full_dataset = full_dataset;
            [GLM_cross_val_shuffle_results.(['file_' file_name])(p).R_square_full_model,...
                GLM_cross_val_shuffle_results.(['file_' file_name])(p).coef_p_val,...
                GLM_cross_val_shuffle_results.(['file_' file_name])(p).R_square_wo_location_model,...                                   % this is redundant - I don't permute the non-location inputs
                GLM_cross_val_shuffle_results.(['file_' file_name])(p).RMSE_full_model,...
                GLM_cross_val_shuffle_results.(['file_' file_name])(p).RMSE_wo_location_model] = ...
                cross_val_prediction(N_blocks, log_vec_perm(:, p), curr_full_dataset, Y, ...
                partial_model_predictor_ind);

        end
    end

end

%%

for i = 1 : length(file_names)
    file_name = file_names{i};
    figure;
    sgtitle(file_name, 'interpreter', 'none')
    subplot(1,3,1)
    histogram([GLM_all_shuffle.(['file_' file_name]).R_square_full_model], 0:0.01:1,...
        'Normalization', 'probability');
    hold on;
    line([GLM_all(f).R_square_full_model,...
        GLM_all(f).R_square_full_model],...
        [0, 1], 'color', 'red')
    title({'R^2 of entire dataset'; 'shuffled location vs. real data'})


    if do_cross_val
        subplot(1,3,2)
        histogram([GLM_cross_val_shuffle_results.(['file_' file_name]).R_square_full_model], 0:0.01:1,...
            'Normalization', 'probability');
        hold on;
        histogram([GLM_cross_val_shuffle_results.(['file_' file_name]).R_square_wo_location_model], 0:0.01:1,...
            'Normalization', 'probability');
        line([GLM_all(f).R_square_full_model,...
            GLM_all(f).R_square_full_model],...
            [0, 1], 'color', 'blue')
        line([GLM_all(f).R_square_wo_location_model  ,...
            GLM_all(f).R_square_wo_location_model  ],...
            [0, 1], 'color', 'red')
        title({'R^2 distribution'; 'of cross validated data'})


        subplot(1,3,3)
        histogram([GLM_cross_val_shuffle_results.(['file_' file_name]).RMSE_full_model], 0:0.5:40,...
            'Normalization', 'probability');
        hold on;
        histogram([GLM_cross_val_shuffle_results.(['file_' file_name]).RMSE_wo_location_model], 0:0.5:40,...
            'Normalization', 'probability');
        line([GLM_all(f).RMSE_full_model,...
            GLM_all(f).RMSE_full_model],...
            [0, 1], 'color', 'blue')
        line([GLM_all(f).RMSE_wo_location_model  ,...
            GLM_all(f).RMSE_wo_location_model  ],...
            [0, 1], 'color', 'red')
        title({'RMSE distribution'; 'of cross validated data'})

    end
end



%% stats and plots - all data, for comparisons A and B:

hist_table = nan(perm_num, length(file_names));
for i = 1 : length(file_names)
    hist_table(:, i) = ...
        [GLM_all_shuffle.(['file_' file_names{i}]).R_square_full_model]';
end

p_val_shuffle_loc = sum(hist_table>[GLM_all.R_square_full_model])./perm_num;

% potential plot for A: likelihood means?

% plot for B:
figure
plot([[GLM_all.R_square_full_model];mean(hist_table)], 'color', [.5 .5 .5])
hold on
plot([[GLM_all.R_square_full_model];mean(hist_table)], '.', 'markersize', 35)
xlim([0 3]);
xticks([1 2])
xticklabels({'real', 'shuffle mean'})
ylabel('r^2')
box off
%% comparison C:
if do_cross_val
hist_table_cross_val_full = nan(perm_num, length(file_names));
hist_table_cross_val_reduced = nan(perm_num, length(file_names));

for i = 1 : length(file_names)
    hist_table_cross_val_full(:, i) = ...
        [GLM_cross_val_shuffle_results.(['file_' file_names{i}]).R_square_full_model]';
    hist_table_cross_val_reduced(:, i) = ...
        [GLM_cross_val_shuffle_results.(['file_' file_names{i}]).R_square_wo_location_model]';
end
% sig_thresh = prctile(hist_table_cross_val, 95);

figure
plot([mean(hist_table_cross_val_full);mean(hist_table_cross_val_reduced)], 'color', [.5 .5 .5])
hold on
plot([mean(hist_table_cross_val_full);mean(hist_table_cross_val_reduced)], '.', 'markersize', 35)
xlim([0 3]);
xticks([1 2])
xticklabels({'full model', 'reduced model'})
ylabel('r^2')
box off

end


%% functions

function [R_square_full_model, R_square_wo_location_model, coefficient_p_val,...
    RMSE_full_model, RMSE_wo_location_model, LLR] = ...
    linear_prediction(full_dataset, Y, partial_model_predictor_ind)

mdl_full = fitlm(zscore(full_dataset), Y);

R_square_full_model = mdl_full.Rsquared.Ordinary;
RMSE_full_model = mdl_full.RMSE;

coefficient_p_val = mdl_full.Coefficients.pValue(2 : end);

mdl_wo_location = fitlm(zscore(full_dataset(:, partial_model_predictor_ind)), Y);

R_square_wo_location_model = mdl_wo_location.Rsquared.Ordinary;
RMSE_wo_location_model = mdl_wo_location.RMSE;

LLR = 2 * (mdl_full.LogLikelihood - mdl_wo_location.LogLikelihood);


end

function log_vec = create_cross_val_log_vec(N_blocks, N_samples)

block_size = floor(1 / N_blocks * N_samples);
block_size_add = mod(N_samples, block_size);
block_size_all = ones(N_blocks, 1) * block_size;
block_size_all(1:block_size_add) = block_size_all(1:block_size_add) + 1;

log_vec = zeros(N_samples, 1);
log_vec(cumsum(block_size_all)+1) = 1;
log_vec = cumsum(log_vec)+1;
log_vec = log_vec(1:end-1);

end

function log_vec = create_cross_val_log_vec_uniform_spread(N_blocks, N_samples)

block_size = floor(1 / N_blocks^2 * N_samples);
block_size_add = rem(N_samples, N_blocks^2);
block_size_all = ones(N_blocks^2, 1) * block_size;
block_size_all(1:block_size_add) = block_size_all(1:block_size_add) + 1;

helper_start = [1; cumsum(block_size_all)];
helper_end = helper_start(2:end) - 1;
helper_end(end) = helper_end(end) + 1;
helper_start = helper_start(1:end-1);

log_vec = zeros(N_samples, 1);

for i = 1 : N_blocks
    for j = 1 : N_blocks
        log_vec(helper_start(i + N_blocks*(j-1)) : helper_end(i + N_blocks*(j-1))) = i;
    end
end

end



function log_vec = ...
    create_log_vec_uniform_with_random_assignment(N_blocks, N_samples,...
    block_size, perms)

% block_size:  X samples in a row! if it's equal to 1, it's the same as the completely random case.

chunks = floor(N_samples/block_size);
if rem(N_samples,block_size) > 0
    chunks = chunks + 1;
end

% not even sampling - some have more samples than others
block_identity = randi(N_blocks, [chunks, perms]);
log_vec = zeros(N_samples, perms);

for p = 1 : perms
    for i = 1 : size(block_identity, 1)
        log_vec((i-1)* block_size + 1 : ...
            min(block_size*i, N_samples), p) = block_identity(i, p);
    end
end
end



function [R_square_full_model, coef_p_val, R_square_wo_loc_model,...
    RMSE_full_model, RMSE_wo_location_model] = ...
    cross_val_prediction(N_blocks, log_vec, full_dataset, Y, ...
    partial_model_predictor_ind)

y_hat = nan(size(log_vec));
y_hat_wo_loc = y_hat;
coef_p_val = nan(size(full_dataset, 2), N_blocks);
for b = 1 : N_blocks
    training_data = ~ismember(log_vec, b);
    testing_data = ~training_data;
    mdl_full = ...
        fitlm(zscore(full_dataset(training_data,:)),...
        Y(training_data));
    y_hat(testing_data) = predict(mdl_full, zscore(full_dataset(testing_data,:)));
    mdl_wo_location = ...
        fitlm(zscore(full_dataset(training_data, partial_model_predictor_ind)),...
        Y(training_data));
    y_hat_wo_loc(testing_data) = ...
        predict(mdl_wo_location,...
        zscore(full_dataset(testing_data,partial_model_predictor_ind)));
    coef_p_val(:, b) = mdl_full.Coefficients.pValue(2:end);
end

R_square_full_model = (corr(y_hat, Y))^2;
R_square_wo_loc_model = (corr(y_hat_wo_loc, Y))^2;

RMSE_full_model = sqrt(mean((y_hat-Y).^2));
RMSE_wo_location_model = sqrt(mean((y_hat_wo_loc-Y).^2));

end