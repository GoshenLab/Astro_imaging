
SAMPLING_RATE = 15.49;
PERM_NUM = 1000;

samples_after_rew = round(SAMPLING_RATE * 7);
discrete_analysis = true;


%% LOAD VARIABLES
% event_mat, df_f, cell_IDs, movement, hardware_ind, lap_vec,
% discrete_location, location_norm


if discrete_analysis
    sig = event_mat;
    concurrent_event_num = sum(sig, 2);
else
    sig = df_F;
    concurrent_event_num = mean(sig, 2, 'omitnan');
end

reward_ind = hardware_ind.valve_ind;

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


p_val_new = sum(R_perm<R_real)./PERM_NUM;

curr_data = after_rew_prop./mean(after_rew_prop(:), 'omitnan');
shadedErrorBar([], mean(curr_data, 'omitnan'),...
    std(curr_data, 'omitnan')./sqrt(sum(~isnan(curr_data))))
hold on;
line([0  size(curr_data, 2)-1], [1 1])

