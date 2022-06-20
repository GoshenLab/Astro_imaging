clear
clc
close all

all_file_names = dir('*_full_sig_ROIs_df_F.signals');

all_quad_file_name = {all_file_names.name};

all_file_names = ...
    cellfun(@(x) x(1:11), all_quad_file_name, 'uniformoutput', false);


for i = 1 : length(all_file_names)
    file_name = all_file_names{i};
    detect_events_both_ways(file_name)
end


function detect_events_both_ways(file_name)

disp(file_name);

if ~isfile([file_name, '_new_event_det_df_F.mat'])
    disp('creating events');
    load([file_name, '_full_sig_ROIs_df_F.signals'], '-mat')

    min_length_thresh = 4;
    plot_me = false;

    R = 1 : size(df_F, 2);
    [final_event_mat, event_hist_bins, event_hist_bins_smoothed] = ...
        get_events_V2(df_F(:, R), plot_me, R);

    % there's a "zero_short_events" input, that's for omitted frames
    [events_above_min_original, event_length, long_event_length] = ...
        get_event_length_above_thresh(final_event_mat, min_length_thresh,...
        true);

    % I added this for movies where there are omitted frames, to
    % replace "no events" between NANs (so there are no "sudden" zeros):
    nan_trim = double(~(sum(isnan(events_above_min_original),2)>0));
    nan_trim_above_min = ...
        get_event_length_above_thresh(nan_trim, min_length_thresh, false);
    events_above_min = events_above_min_original;
    events_above_min(isnan(nan_trim_above_min), :) = nan;

    save([file_name, '_new_event_det_df_F.mat'],...
        'final_event_mat', 'event_hist_bins',...
        'event_hist_bins_smoothed', 'events_above_min', 'event_length', ...
        'long_event_length', 'min_length_thresh');

else
    disp('events already created');

end


end
