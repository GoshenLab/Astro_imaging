function [events_above_min, event_length, long_event_length] = ...
    get_event_length_above_thresh(events, min_length_thresh, zero_short_events)

% inputs:
%   events:             event matrix...
%   min_length_thresh:  in samples, not time
%   zero_short_events:  turn below threshold to zero (true) or NAN (false)-
%                       this is to enable running this function on vectors
%                       with nans

% outputs:
%   events_above_min:   logical mat, only long enough events.
%   event_length:       all events length (before omitting)
%   long_event_length:  only events above length "min_length_thresh"

SAMPLING_RATE = 15.49;
event_length = cell(size(events, 2), 1);
long_event_length = event_length;
events_above_min = events;

for s = 1 : size(events, 2)
    % if there is no 0 at the begining/end, you will miss the event
    % (logically - it might be trimmed)
    f = find(diff([false; events(:, s) > 0; false]) ~= 0);
    event_length{s} =...
        [f(2 : 2 : end) - f(1 : 2 : (end - 1))]./SAMPLING_RATE;
    long_event_length{s} = ...
        event_length{s}(event_length{s} > (min_length_thresh/SAMPLING_RATE));
    short_events_ind = ...
        find(event_length{s} <= (min_length_thresh/SAMPLING_RATE));
    if ~isempty(short_events_ind)
        startush = f(1 : 2 : (end - 1));
        endush = f(2 : 2 : end) - 1;
        for i = short_events_ind'
            if zero_short_events
                events_above_min(startush(i) : endush(i), s) = 0;
            else
                events_above_min(startush(i) : endush(i), s) = nan;
            end
        end
    end
end


%%

end