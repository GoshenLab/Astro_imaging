
quad_or_vel = 1; % 1: location, 2: velocity

conf_bin = 0.1;

train_prop = 0.8;
sim_num = 1000;

%% binning for histogram
bin_size = 0.02;
bin_vec = [0 : bin_size : 1];
H = nan(4, length(bin_vec) - 1, sim_num);

%% load variables:
% df_F, event_mat, hardware_ind, movement, location_norm, location_discrete,
% lap_vec, velocity, ROI_IDs

%% prepare signal:
curr_sig_original = event_mat;  % or: dF/F
curr_sig_original = curr_sig_original(movement == 1, :);

if quad_or_vel == 1
    curr_data = location_norm(movement == 1);
    curr_data(curr_data < 0) = 0;
    curr_data(curr_data > 1) = 1;
else
    curr_data = velocity(movement == 1);
    curr_data = ...
        normalize(curr_data', 'range');
end

lap_vec = lap_vec(movement == 1);

test_length = length(curr_data) - ...
    round(train_prop * length(curr_data));

error_real_pooled = nan(test_length, sim_num);
error_perm_pooled = error_real_pooled;

conf_mat_real = nan(1/conf_bin, 1/conf_bin, sim_num);
conf_mat_perm = conf_mat_real;

perm_delta = randperm(length(curr_data), sim_num);



for i = 1 : sim_num
    if mod(i, 100) == 0
        disp(i)
    end
    [train_ind, test_ind] = ...
        get_train_test_inds(train_prop, length(curr_data));

    y_train = curr_data(train_ind)';
    y_test = curr_data(test_ind)';
    y_test_dis = categorical(discretize(y_test, [0 : conf_bin : 1]));


    curr_sig = curr_sig_original;

    [~, y_hat_test] = ...
        do_regression(train_ind, test_ind, curr_sig, y_train);

    %% trimming (relevant for both velocity and location):

    y_hat_test = min(y_hat_test, 1);
    y_hat_test = max(y_hat_test, 0);
    y_hat_test_dis = ...
        categorical(discretize(y_hat_test, [0 : conf_bin : 1]));

    error_real = abs(y_hat_test - y_test);

    C = confusionmat(y_test_dis, y_hat_test_dis);

    conf_mat_real(:, :, i) = C;

    %% permutation
    sig_perm = ...
        circshift(curr_sig_original, perm_delta(i));

    [~, y_hat_test_perm] = ...
        do_regression(train_ind, test_ind, sig_perm, y_train);

    y_hat_test_perm = min(y_hat_test_perm, 1);
    y_hat_test_perm = max(y_hat_test_perm, 0);

    y_hat_test_perm_dis = ...
        categorical(discretize(y_hat_test_perm, [0 : conf_bin : 1]));

    C = confusionmat(y_test_dis, y_hat_test_perm_dis);

    if length(conf_mat_perm(:, :, i)) == length(C)
        conf_mat_perm(:, :, i) = C;
    end

    error_perm = abs(y_hat_test_perm - y_test);

    h_cdf_real = histogram(error_real, bin_vec,...
        'Normalization', 'cdf',...
        'Visible', 'off');
    h_cdf_real = h_cdf_real.Values;

    h_cdf_perm = histogram(error_perm, bin_vec,...
        'Normalization', 'cdf',...
        'Visible', 'off');
    h_cdf_perm = h_cdf_perm.Values;

    H(1, :, i) = h_prob_real;
    H(2, :, i) = h_prob_perm;
    H(3, :, i) = h_cdf_real;
    H(4, :, i) = h_cdf_perm;

    H_mean = mean(H, 3);
    H_SE = std(H, [], 3)./sqrt(sim_num);
    subplot(1,2,1)
    shadedErrorBar(bin_vec(1:end-1) + bin_size/2,...
        H_mean(1, :), H_SE(1, :), 'lineprops', 'b');
    hold on;
    shadedErrorBar(bin_vec(1:end-1) + bin_size/2,...
        H_mean(2, :), H_SE(2, :), 'lineprops', 'r');
    xlabel('error size')
    ylabel('probability')

    subplot(1,2,2)
    line([0 1], [0 1], 'LineStyle', '--', 'Color', 'k')
    hold on;
    shadedErrorBar(bin_vec(1:end-1) + bin_size/2,...
        H_mean(3, :), H_SE(3, :), 'lineprops', 'b');
    shadedErrorBar(bin_vec(1:end-1) + bin_size/2,...
        H_mean(4, :), H_SE(4, :), 'lineprops', 'r');
    xlabel('error size')
    ylabel('cumulative probability')


    error_real_pooled(:, i) = error_real;
    error_perm_pooled(:, i) = error_perm;

end


%% stats:
p_val = ...
    1 - (...
    sum(mean(error_real_pooled) <= mean(error_perm_pooled))./...
    sim_num);

if quad_or_vel == 1
    norm_2_belt = 170;
else
    norm_2_belt = 1;
end


mean_real = mean(error_real_pooled * norm_2_belt,'omitnan');
mean_perm = mean(error_perm_pooled * norm_2_belt, 'omitnan');

overall_mean_real = mean(mean_real);
overall_mean_perm = mean(mean_perm);

SE_real = ...
    (std(mean_real, 'omitnan')./sqrt(sum(~isnan(mean_real))) * norm_2_belt);
SE_perm = ...
    (std(mean_perm, 'omitnan')./sqrt(sum(~isnan(mean_perm))) * norm_2_belt);


%% confusion mat:
h_f = figure();
h_ax1 = subplot(2,1,1);
imagesc(mean(conf_mat_real./sum(conf_mat_real, 2), 3, 'omitnan'))
if quad_or_vel == 1
    ylabel('real location')
    xlabel('predicted location')
else
    ylabel('real velocity')
    xlabel('predicted velocity')
end
title('real data')
colorbar()
axis equal
h_ax1.TickDir = 'out';
h_ax1.XLim = [0.5, 1/conf_bin + 0.5];

h_ax2 = subplot(2,1,2);
imagesc(mean(conf_mat_perm./sum(conf_mat_perm, 2), 3, 'omitnan'))
if quad_or_vel == 1
    ylabel('real location')
    xlabel('predicted location')
else
    ylabel('real velocity')
    xlabel('predicted velocity')
end
title('shuffled data')
colorbar()
axis equal
h_ax2.TickDir = 'out';
h_ax2.XLim = [0.5, 1/conf_bin + 0.5];

max_val =...
    max([max(h_ax1.Children.CData(:)), max(h_ax2.Children.CData(:))]);

h_ax1.CLim(2) = max_val;
h_ax2.CLim(2) = max_val;
h_ax1.CLim(1) = 0;
h_ax2.CLim(1) = 0;





%% functions:

function [y_hat_train, y_hat_test] = ...
    do_regression(train_ind, test_ind, sig, y_train)

train_sig = [ones(length(train_ind), 1), sig(train_ind, :)];
test_sig = sig(test_ind, :);

w_hat = train_sig \ y_train;
y_hat_train = (train_sig * w_hat);
y_hat_test = [w_hat(1) + test_sig * w_hat(2 : end)];
end



