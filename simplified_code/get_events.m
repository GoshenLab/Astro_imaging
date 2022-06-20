function [final_event_mat, event_hist_bins, event_hist_bins_smoothed] = ...
    get_events_V2(deltaF_dombeck_sig_df_f, plot_me, ROI_names)


% good frames = legacy. Now I use a constant number of bins, independent of
% signal length
% R: ROI names


SAMPLING_RATE = 15.49;
ROWS = 5;
COL = 2;

if plot_me
    c = 1;
end

% get events based on the original df/F traces
event_hist_bins = 30;
[event_mat, ~] = ...
    find_event_threshold(deltaF_dombeck_sig_df_f, event_hist_bins);

% smooth the df/F traces, and get the threshold! (not the event mat)
df_F_mat_smooth = smoothdata(deltaF_dombeck_sig_df_f, 'movmedian', 61);
event_hist_bins_smoothed = 75;
df_F_mat_smooth(isnan(deltaF_dombeck_sig_df_f)) = nan;
[event_mat_smoothed, threshold_smoothed] = ...
    find_event_threshold(df_F_mat_smooth, event_hist_bins_smoothed);
event_mat_smoothed(isnan(deltaF_dombeck_sig_df_f)) = nan;
% threshold_mean = mean([threshold; threshold_smoothed]);

% I canceled the next lines because you can calculate this in the function
% event_mat_smoothed = double(df_F_mat_smooth > threshold_smoothed);
% event_mat_smoothed(isnan(deltaF_dombeck_sig_df_f)) = nan;
 
% create mat - 1 is when both detection methods find an event
event_mat_both = double((event_mat_smoothed + event_mat) == 2);
event_mat_both(isnan(event_mat)) = nan;

% this is only relevant when there are no NANs (I turn them to zeros); I
% look for differences, and search for positive events, so this is not
% calculated. I use the original matrix, so the final product still has NAN
test_event_mat_smoothed = event_mat_smoothed;
test_event_mat_smoothed((isnan(test_event_mat_smoothed))) = 0;

event_mat_both_copy = event_mat_both;
event_mat_both_copy(isnan(event_mat_both)) = 0;

test_event_mat_smoothed = [zeros(1, size(test_event_mat_smoothed, 2));...
    test_event_mat_smoothed;...
    zeros(1, size(test_event_mat_smoothed, 2))];

final_event_mat = event_mat_smoothed;

[event_start, cell_id] = find(diff(test_event_mat_smoothed) == 1);
[event_end, ~] = find(diff(test_event_mat_smoothed) == -1);
event_end = event_end - 1;

for i = unique(cell_id)'
    curr_start_ind = event_start(cell_id == i);
    curr_end_ind = event_end(cell_id == i);
    for j = 1 : length(curr_start_ind)
        % how many bins are double using both event detection methods):
        test = ...
            sum(event_mat_both_copy(curr_start_ind(j) : curr_end_ind(j), i) == 1); 
        % if there isn't a "simple" event detection within smoothed, cancel
        % the event
        if test == 0
            final_event_mat(curr_start_ind(j) : curr_end_ind(j), i) = 0;
        end
    end
end

%%
color_map = lines(3);
if plot_me
    figure;
    x = [1 : size(deltaF_dombeck_sig_df_f, 1)]./ (SAMPLING_RATE * 60);
    for i = 1 : size(deltaF_dombeck_sig_df_f, 2)
        subplot(ROWS, COL, c)
        
        %% threshold based on regular/smoothed data:
        plot(x, deltaF_dombeck_sig_df_f(:, i));
%         condA = deltaF_dombeck_sig_df_f(:, i) > threshold(i);
%         condB = df_F_mat_smooth(:, i) > threshold_smoothed(i);
%         hold on;
%         plot(x(condB),...
%             deltaF_dombeck_sig_df_f(condB, i), 'r.', 'MarkerSize', 4)
%         plot(x(condA),...
%             deltaF_dombeck_sig_df_f(condA, i), 'y.', 'MarkerSize', 2);
%         plot(x(condA & condB),...
%             deltaF_dombeck_sig_df_f(condA & condB, i), 'k.');

hold on;
temp = find(final_event_mat(:, i));
plot(x(temp), deltaF_dombeck_sig_df_f(temp, i), 'k.');
        %% threshold based on combined:
        
        %         plot(x, deltaF_dombeck_sig_df_f(:,i));
        %         condA = deltaF_dombeck_sig_df_f(:,i) > threshold(i);
        %         condB = df_F_mat_smooth(:,i) > threshold_mean(i);
        %         hold on;
        %         plot(x(condB),...
        %             deltaF_dombeck_sig_df_f(condB,i), '.', 'MarkerSize', 4)
        %         plot(x(condA),...
        %             deltaF_dombeck_sig_df_f(condA,i), '.');
        %%
        
        
        
        
        title(['cell num: ' num2str(ROI_names(i))])
        xlabel('time (min)')
        c = c + 1;
        if c > ROWS * COL
            figure;
            c = 1;
        end
    end
end





end



function [event_mat, threshold] = ...
    find_event_threshold(sig, event_hist_bins)


event_mat = nan(size(sig));

threshold = zeros(1, size(sig, 2));

for i = 1 : size(sig, 2)
    % since the smoothed version fills in Nans, omit them based on the
    % original signal:
    included_frames = ~isnan(sig(:, i));
    
    a_h = histogram(sig(included_frames,i), event_hist_bins,...
        'Visible', 'off');
    [~, t] = max(a_h.Values);
    % get the mode value od the histogram
    mode_val = (a_h.BinEdges(2) - a_h.BinEdges(1))/2 + a_h.BinEdges(t);
    % get the minimal value, based on the the mean of the lowest 100
    % samples
    min_val = mean(mink(sig(included_frames,i), 100));
    % the definition of the threshold
    threshold(i) = (2 * mode_val - min_val);
    
    event_mat(included_frames, i) =...
        sig(included_frames, i) > threshold(i);
    
    
end
end