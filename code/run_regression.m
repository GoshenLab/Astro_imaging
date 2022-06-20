%%
quad_or_vel = 1; % 1: location, 2: velocity
conf_bin = 0.1;

SAMPLING_RATE = 15.49;
train_prop = 0.8;
sim_num = 10;

save_here = ['D:\astro_imaging\Nature_code\all_code_data\',...
    'decoder_location_results_',...
    'train_prop_', num2str(train_prop), '\'];

if ~isfolder(save_here)
    error('no save here folder')
end

%%
data_folder = 'D:\astro_imaging\Nature_code\all_code_data\';

%% familiar (location or velocity)
file_names = {...
    '8C3_100_001',...
    '9H3_003_000',...
    '9I3_004_000',...
    '9P2_000_000',...
    '9N1_000_000',...
    '9Q4_000_000',...
    };


%% novel:
% file_names = {...
%     '8C3_000_000';
%     '9B2_000_002';
%     '9P2_001_000';
%     '9N1_001_000';
%     '9O3_001_000';
%     '9Q4_001_000';
%     '9S3_001_000';
%     };


H = cell(size(file_names));
error_real_pooled = H;
error_perm_pooled = H;
conf_mat_perm = H;
conf_mat_real = H;

error_struct_cell = H;

for i = 1 : length(file_names)
    disp(file_names{i})
    if ~isfile([save_here, file_names{i}, '_error_pooled_', num2str(sim_num), '.mat'])
        disp('running decoder')
        [temp_error_real_pooled, temp_error_perm_pooled,...
            temp_conf_mat_real, temp_conf_mat_perm,...
            temp_error_struct_cell] = ...
            regression_clean(file_names{i},...
            train_prop, sim_num,...
            quad_or_vel, conf_bin...
            );

        save([save_here, file_names{i}, '_error_pooled_', num2str(sim_num), '.mat'],...
            'temp_error_real_pooled',...
            'temp_error_perm_pooled',...
            'temp_error_struct_cell',...
            'temp_conf_mat_real',...
            'temp_conf_mat_perm')

    else
        disp('loading results')
        load([save_here, file_names{i}, '_error_pooled_', num2str(sim_num), '.mat'],...
            'temp_error_real_pooled',...
            'temp_error_perm_pooled',...
            'temp_error_struct_cell',...
            'temp_conf_mat_real',...
            'temp_conf_mat_perm'...
            )
    end
    error_real_pooled{i} = temp_error_real_pooled;
    error_perm_pooled{i} = temp_error_perm_pooled;
    conf_mat_real{i} = temp_conf_mat_real;
    conf_mat_perm{i} = temp_conf_mat_perm;
    error_struct_cell{i} = temp_error_struct_cell;
end



%%
BELT_LENGTH = 170;

bin_size = 0.02;

if quad_or_vel == 1
    bin_vec = [0 : bin_size : 1] * BELT_LENGTH;
    error_real_pooled = ...
        cellfun(@(x) x * BELT_LENGTH, error_real_pooled,...
        'UniformOutput', false);
    error_perm_pooled = ...
        cellfun(@(x) x * BELT_LENGTH, error_perm_pooled,...
        'UniformOutput', false);

else
    bin_vec = [0 : bin_size : 1];
end

%% histograms for pooled data!
h_cdf_real = nan(length(bin_vec) - 1, length(file_names));
h_cdf_perm = h_cdf_real;
h_cdf_real_mode = h_cdf_real;

h_f = figure();
COLORS = lines(length(file_names));

% hold on;
for i = 1 : length(file_names)
    subplot(2,5,i)

    temp = histogram(error_real_pooled{i}, bin_vec,...
        'Normalization', 'cdf',...
        'Visible', 'off');
    h_cdf_real(:, i) = temp.Values;
    plot(bin_vec(2 : end), h_cdf_real(:, i), 'Color', COLORS(i, :))
    hold on
    temp = histogram(error_perm_pooled{i}, bin_vec,...
        'Normalization', 'cdf',...
        'Visible', 'off');
    h_cdf_perm(:, i) = temp.Values;
    plot(bin_vec(2 : end), h_cdf_perm(:, i), '--', 'Color', COLORS(i, :))

    xlim([0 max(bin_vec)])
    if quad_or_vel == 1
        xlabel('error size (cm)')
        xticks([0 85 170])

    else
        xlabel('error size (relative velocity)')
        xticks([0 0.5 1])
    end
    ylabel('cumulative probability')
    title(file_names{i}, 'interpreter', 'none')
end


subplot(2,5,10)
shadedErrorBar(bin_vec(2 : end),...
    mean(h_cdf_real, 2, 'omitnan'),...
    std(h_cdf_real, [], 2, 'omitnan')./sqrt(sum(~isnan(h_cdf_real), 2)),...
    'lineprops', 'b');
hold on;
shadedErrorBar(bin_vec(2 : end),...
    mean(h_cdf_perm, 2, 'omitnan'),...
    std(h_cdf_perm, [], 2, 'omitnan')./sqrt(sum(~isnan(h_cdf_real), 2)),...
    'lineprops', 'r');
if quad_or_vel == 1
    xlabel('error size (cm)')
    xlim([0 max(bin_vec)])
    xticks([0 85 170])
else
    xlabel('error size (relative velocity)')
end
ylabel('cumulative probability')


%% stats

mean_real = cellfun(@mean, error_real_pooled, 'UniformOutput', false);
mean_real = cell2mat(cellfun(@mean, mean_real, 'UniformOutput', false));

mean_perm = cellfun(@mean, error_perm_pooled, 'UniformOutput', false);
mean_perm = cell2mat(cellfun(@mean, mean_perm, 'UniformOutput', false));

disp('mean real error')
disp(mean(mean_real, 'omitnan') );
disp('SE real error')
disp(std(mean_real, 'omitnan')./sqrt(sum(~isnan(mean_real))) );

disp('mean perm error')
disp(mean(mean_perm, 'omitnan') );
disp('SE real error')
disp(std(mean_perm, 'omitnan')./sqrt(sum(~isnan(mean_perm))) );

%% 
p_val = nan(length(file_names), 1);
for i = 1 : length(file_names)
    p_val(i) = ...
        1 - (...
        sum(mean(error_real_pooled{i}) <= mean(error_perm_pooled{i}))./...
        sim_num);
end

%% confusion matrices
% mean probability, averaging across permutations

for i = 1 : length(file_names)
    h_f(i) = figure();
    sgtitle(file_names{i}, 'interpreter', 'none')
    h_ax1 = subplot(2,1,1);
    imagesc(mean(conf_mat_real{i}./sum(conf_mat_real{i}, 2), 3, 'omitnan'))
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
    imagesc(mean(conf_mat_perm{i}./sum(conf_mat_perm{i}, 2), 3, 'omitnan'))
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

end





%% functions:

function [error_real_pooled, error_perm_pooled,...
    conf_mat_real, conf_mat_perm,...
    error_struct] = ...
    regression_clean(file_name,...
    train_prop, sim_num,...
    quad_or_vel, conf_bin)

mouse_name = file_name(1:3);

bad_test_set = 0;

%% errors:
L1_error_real_trim = nan(sim_num, 1);
L1_error_perm_trim = L1_error_real_trim;

bin_size = 0.02;
bin_vec = [0 : bin_size : 1];
H = nan(4, length(bin_vec) - 1, sim_num);


[df_F, events_above_min, movement,...
    quad_data_norm, lap_vec_custom,...
    velocity_by_quad_data_smooth] = ...
    get_variables_regression(file_name);

if isfile([file_name, '_corr_ind.mat'])
    corr_ind = sum(isnan(events_above_min), 2) == 0;
    events_above_min = events_above_min(corr_ind, :);
    df_F = df_F(corr_ind, :);
    movement = movement(corr_ind);
    quad_data_norm = quad_data_norm(corr_ind);
    lap_vec_custom = lap_vec_custom(corr_ind);
    velocity_by_quad_data_smooth = velocity_by_quad_data_smooth(corr_ind);
else
    corr_ind = ones(size(movement));
end


%% final datasets:

curr_sig_original = events_above_min;
curr_sig_original = curr_sig_original(movement == 1, :);
if quad_or_vel == 1
    curr_data = quad_data_norm(movement == 1);
    curr_data(curr_data < 0) = 0;
    curr_data(curr_data > 1) = 1;
else
    curr_data = velocity_by_quad_data_smooth(movement == 1);
    curr_data = ...
        normalize(curr_data', 'range');
end
lap_vec_custom = lap_vec_custom(movement == 1);

bad_ind = ...
    (isnan(curr_data))' | (sum(isnan(curr_sig_original), 2)>0);

if sum(bad_ind) > 0
    disp('something is wrong - there are NANs!')
    curr_sig_original(bad_ind, :) = [];
    curr_data(bad_ind) = [];
    lap_vec_custom(bad_ind) = [];
end

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

    [y_hat_train, y_hat_test] = ...
        do_regression(train_ind, test_ind, curr_sig, y_train);

    %% trimming (relevant for both velocity and location):

    y_hat_test = min(y_hat_test, 1);
    y_hat_test = max(y_hat_test, 0);
    y_hat_test_dis = ...
        categorical(discretize(y_hat_test, [0 : conf_bin : 1]));

    error_real = abs(y_hat_test - y_test);
    L1_error_real_trim(i) = mean(error_real);

    C = confusionmat(y_test_dis, y_hat_test_dis);
    if length(conf_mat_real(:, :, i)) ~= length(C)
        disp('not enough samples')
        bad_test_set = bad_test_set + 1;
    else
        conf_mat_real(:, :, i) = C;
    end


    %% permutation!
    sig_perm = ...
        circshift(curr_sig_original, perm_delta(i));

    [~, y_hat_test_perm] = ...
        do_regression(train_ind, test_ind, sig_perm, y_train);


    %% trimming (relevant for both velocity and location):

    y_hat_test_perm = min(y_hat_test_perm, 1);
    y_hat_test_perm = max(y_hat_test_perm, 0);

    y_hat_test_perm_dis = ...
        categorical(discretize(y_hat_test_perm, [0 : conf_bin : 1]));

    C = confusionmat(y_test_dis, y_hat_test_perm_dis);

    if length(conf_mat_perm(:, :, i)) == length(C)
        conf_mat_perm(:, :, i) = C;
    end

    error_perm = abs(y_hat_test_perm - y_test);
    
    L1_error_perm_trim(i) = mean(error_perm);


    h_prob_real = histogram(error_real, bin_vec,...
        'Normalization', 'probability',...
        'Visible', 'off');
    h_prob_real = h_prob_real.Values;
    h_prob_perm = histogram(error_perm, bin_vec,...
        'Normalization', 'probability',...
        'Visible', 'off');
    h_prob_perm = h_prob_perm.Values;

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

    error_real_pooled(:, i) = error_real;
    error_perm_pooled(:, i) = error_perm;

end

disp(mouse_name)
disp('bad test set number')
disp(bad_test_set)

%%
figure();
H_mean = mean(H, 3);
H_SE = std(H, [], 3)./sqrt(sim_num);
sgtitle(file_name, 'interpreter', 'none')
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



%% store errors in struct:
error_struct = struct();

error_struct.L1_error_real_cyc = L1_error_real_trim;
error_struct.L1_error_perm_cyc = L1_error_perm_trim;

error_struct.histogram_data = H;

end


function [y_hat_train, y_hat_test] = ...
    do_regression(train_ind, test_ind, sig, y_train)

train_sig = [ones(length(train_ind), 1), sig(train_ind, :)];
test_sig = sig(test_ind, :);

w_hat = train_sig \ y_train;
y_hat_train = (train_sig * w_hat);
y_hat_test = [w_hat(1) + test_sig * w_hat(2 : end)];
end

