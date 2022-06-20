
min_length_thresh = 4;
plot_me = false;

R = 1 : size(df_F, 2);
[final_event_mat, event_hist_bins, event_hist_bins_smoothed] = ...
    get_events(df_F(:, R), plot_me, R);

zero_short_events = true;
% there's a "zero_short_events" input, that's for omitted frames
[events_above_min_original, event_length, long_event_length] = ...
    get_event_length_above_thresh(final_event_mat, min_length_thresh,...
    zero_short_events);

% I added this for movies where there are omitted frames, to
% replace "no events" between NANs (so there are no "sudden" zeros):
nan_trim = double(~(sum(isnan(events_above_min_original),2)>0));
nan_trim_above_min = ...
    get_event_length_above_thresh(nan_trim, min_length_thresh, false);
events_above_min = events_above_min_original;
events_above_min(isnan(nan_trim_above_min), :) = nan;
