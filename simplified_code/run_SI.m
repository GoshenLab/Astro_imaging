
sim_num = 1000;
discrete_analysis = true; % or false for df/f

%% LOAD VARIABLES
% signals (df/f or event_mat, ROI_IDs, lap_vec, discrete_location,
% movement, velocity, hardware_ind

%%
good_frames = ...
    (sum(isnan(sig), 2) == 0) & (movement == 1) &...
    ~isnan(discrete_location');

sig = sig(good_frames, ROI_IDs);

discrete_location = discrete_location(good_frames);

time_spent_in_bin = ...
    splitapply(@length, discrete_location, discrete_location);

time_spent_in_bin = time_spent_in_bin(2:end);

total_time = sum(time_spent_in_bin);

p_bin = [time_spent_in_bin./total_time];

event_prob_in_bin = ...
    splitapply(@mean, sig, discrete_location');
event_prob_in_bin = event_prob_in_bin(2:end, :);

mean_event_prob = p_bin * event_prob_in_bin;

norm_event_prob = event_prob_in_bin./mean_event_prob;
norm_event_prob(isequal(norm_event_prob, inf)) = nan;
SI_real = ...
    sum(p_bin' .* norm_event_prob ...
    .* log2(norm_event_prob), 'omitnan');

perm_delta = randperm(length(sig), sim_num);
SI_perm = nan(sim_num, size(sig, 2));
for sim = 1 : sim_num
    if mod(sim, 100) == 0
        disp(['sim num: '  num2str(sim)])
    end
    curr_sig = circshift(sig, perm_delta(sim));
    event_prob_in_bin = ...
        splitapply(@mean, curr_sig, discrete_location');
    event_prob_in_bin = event_prob_in_bin(2: end, :);
    mean_event_prob = p_bin * event_prob_in_bin;
    norm_event_prob_perm = event_prob_in_bin./mean_event_prob;
    norm_event_prob_perm(isequal(norm_event_prob_perm, inf)) = nan;
    SI_perm(sim, :) = ...
        sum(p_bin' .* norm_event_prob_perm ...
        .* log2(norm_event_prob_perm), 'omitnan');
end

sig_SI = SI_real > prctile(SI_perm, 95);


