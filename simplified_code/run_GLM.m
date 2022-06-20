
partial_model_predictor_ind = [2:5]; % the columns to include (1 = location)
perm_num = 1000;

%%
% LOAD VARIABLES: lap_vec, discrete_location, velocity, event_mat,
% hardware_ind, movement, ROIs
%%
total_lap_number = max(lap_vec);
elapsed_time = nan(size(lap_vec)); % in samples (starting at zero)
for i = 1 : total_lap_number
    curr_ind = find(lap_vec == i);
    elapsed_time(curr_ind : end) = ...
        0 : length(curr_ind : length(elapsed_time) - 1);
end

%%
win_length = 30;
integration_filter = ones(win_length, 1);
velocity_integral = ...
    conv2(velocity, integration_filter, 'same');

%%
licking = hardware_ind.lick_ind;
licking_vec = zeros(size(lap_vec_custom));
licking_vec(licking) = 1;

%%
concurrent_events = sum(event_mat, 2);
good_sig = ~isnan(concurrent_events);   % for bad registration removal

%%

X = [...
    discrete_location',...
    elapsed_time',...
    velocity,...
    velocity_integral,...
    licking_vec'...
    ];

include_frames = ...
    find((sum(isnan(X),2) == 0) & (movement == 1) & good_sig);


full_dataset = X(include_frames, :);
Y = concurrent_events(include_frames);

N_samples = length(full_dataset);

%% testing on all the dataset, wo cross val:

[R_square_full_model,...
    R_square_wo_location_model,...
    coefficient_p_val,...
    RMSE_full_model,...
    RMSE_wo_location_model,...
    LLR] = ...
    linear_prediction(full_dataset, Y, partial_model_predictor_ind);

% shuffle the location vector:
perm_ind = randi(N_samples, [perm_num, 1]);

R_square_full_model_shuffle = nan(perm_num, 1);
coefficient_p_val_shuffle = R_square_full_model_shuffle;
RMSE_full_model_shuffle = R_square_full_model_shuffle;
RMSE_wo_location_model_shuffle = R_square_full_model_shuffle;

for p = 1 : perm_num
    if mod(p, 100) == 0
        disp(p)
    end
    curr_full_dataset = full_dataset;
    curr_full_dataset(:, 1) = ...            % LOCATION is in the first column!
        circshift(full_dataset(:, 1), perm_ind(p));
    [R_square_full_model_shuffle,...
        ~,...                                % I don't permute any input, so it's redundant
        coefficient_p_val_shuffle,...
        RMSE_full_model_shuffle,...
        RMSE_wo_location_model_shuffle] = ...
        linear_prediction(curr_full_dataset, Y, partial_model_predictor_ind);
end


%% cross-validation
N_blocks = 10;
block_size = 31; % 2 sec
log_vec_perm = ...
    create_log_vec_uniform_with_random_assignment(...
    N_blocks, N_samples, block_size, perm_num);

R_square_full_model_Xval = nan(perm_num, 1);
coef_p_val_Xval = R_square_full_model_Xval;
R_square_wo_location_model_Xval = R_square_full_model_Xval;
RMSE_full_model_Xval = R_square_full_model_Xval;
RMSE_wo_location_model_Xval = R_square_full_model_Xval;

for p = 1 : perm_num
    if mod(p, 100) == 0
        disp(p)
    end

    curr_full_dataset = full_dataset;
    [R_square_full_model_Xval(p),...
        coef_p_val_Xval(p),...
        R_square_wo_location_model_Xval(p),...
        RMSE_full_model_Xval(p),...
        RMSE_wo_location_model_Xval(p)] = ...
        cross_val_prediction(N_blocks, log_vec_perm(:, p),...
        curr_full_dataset, Y, ...
        partial_model_predictor_ind);

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


function log_vec = ...
    create_log_vec_uniform_with_random_assignment(N_blocks, N_samples,...
    block_size, perms)

% block_size:  X samples in a row! if it's equal to 1, it's the same as the completely random case.

chunks = floor(N_samples/block_size);
if rem(N_samples,block_size) > 0
    chunks = chunks + 1;
end

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