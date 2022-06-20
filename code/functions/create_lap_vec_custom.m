function [lap_vec_custom, discrete_quad_data, quad_data_norm] = ...
    create_lap_vec_custom(quad_data, bin_number, hardware_ind,...
    remove_us_vec)

% INPUTS:
% quad_data:    the raw quadrature file, without modulu function! use after
%               trimming the first frames.
% hardware_ind: the 'moved' hardware, from the saved variable (after trim)
% remove_us_vec: the trimmed variable!

% OUTPUTS:
% lap_vec_custom        division into laps, using the bin number
% discrete_quad_data    discrete values (1:bin_number) - if the mouse moved
%                       backwards, you will get the modulu values
% quad_data_norm        normalized values (0:1) - if the mouse moved
%                       backwards, you'll get negative values

% in this function, I calculate the relevant location within lap such that
% the number of ticks per lap is not constant! a lap is defined as:
% (rew_ind + 1) : rew_ind (you omit incomplete laps):


quad_data = double(quad_data);
quad_data(find(remove_us_vec)) = nan;

valve_ind = hardware_ind.valve_ind;

quad_data_norm = nan(size(quad_data));
discrete_quad_data = quad_data_norm;
lap_vec_custom = zeros(size(quad_data_norm));

for i = 1 : (length(valve_ind) - 1)
    curr_ind = valve_ind(i) : valve_ind(i + 1);
    curr_quad = quad_data(curr_ind);
    if sum(isnan(curr_quad)) > 0 && ~isnan(curr_quad(1))
        curr_quad(1) = nan;
        min_quad_ind = find(~isnan(curr_quad), 1);
    elseif sum(isnan(curr_quad)) == 1 && isnan(curr_quad(1))
        min_quad_ind = find(~isnan(curr_quad), 1);
    else
        min_quad_ind = 1;
    end
    
    % added condition for completely removed laps
    if sum(isnan(curr_quad)) < length(curr_quad)
        
        quad_data_norm(curr_ind) = ...
            (curr_quad - curr_quad(min_quad_ind))./...
            (curr_quad(end) - curr_quad(min_quad_ind));
        % subtract the first quad data non-nan value (not ness. the minimal
        % value!!)
        discrete_quad_data(curr_ind) = ...
            discretize(mod(quad_data_norm(curr_ind), max(quad_data_norm(curr_ind))), ...
            linspace(0, 1, bin_number + 1));
        
    
    else
        disp(i)
    end

    lap_vec_custom(valve_ind(i)) = 1;
end
lap_vec_custom(valve_ind(end)) = 1;

% if there's a bad lap in the middle of a run, it will be counted but
% empty:
lap_vec_custom = cumsum(lap_vec_custom, 'omitnan');
lap_vec_custom(isnan(quad_data_norm)) = nan;

%% visualize lap devision:
figure;
subplot(1,2,1)
scatter(1 : length(quad_data_norm), quad_data_norm, 5,...
    discrete_quad_data, 'filled');
colormap(gca, 'hsv')
hold on;
plot(valve_ind, quad_data_norm(valve_ind), 'kX');
IR_ind = setdiff(hardware_ind.IR_ind, find(remove_us_vec==1));
plot(IR_ind, quad_data_norm(IR_ind), 'k^');
legend('path', 'valve', 'IR')
subplot(1,2,2)
scatter(1 : length(quad_data_norm), quad_data_norm, 5,...
    lap_vec_custom, 'filled');
colormap(gca, 'colorcube')

end