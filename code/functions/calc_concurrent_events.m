function [concurrent_event_num, template_mat] = ...
    calc_concurrent_events(movement, sig, lap_vec,...
    discrete_quad_data_shift, discrete_analysis)


movement = logical(movement);
if discrete_analysis
    concurrent_event_num = sum(sig, 2);
else
    concurrent_event_num = mean(sig, 2);
end

template_mat = number_of_sim_events_in_bin_per_lap(lap_vec(movement),...
    discrete_quad_data_shift(movement), concurrent_event_num(movement), 1);

end


