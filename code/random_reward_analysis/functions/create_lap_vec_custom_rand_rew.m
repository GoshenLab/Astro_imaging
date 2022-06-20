function [lap_vec_custom, discrete_quad_data, quad_data_norm, remove_us_vec] = ...
    create_lap_vec_custom_rand_rew(quad_data, bin_number, hardware_ind,...
    remove_us_vec, IR_shift)

% INPUTS:
% quad_data:    the raw quadrature file, without modulu function! use after
%               trimming the first frames. 
% hardware_ind: the 'moved' hardware, from the saved variable
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


%IR_ind = setdiff(hardware_ind.IR_ind, find(remove_us_vec==1));
IR_ind = hardware_ind.IR_ind;
valve_ind = hardware_ind.valve_ind;
quad_data_norm = nan(size(quad_data));
discrete_quad_data = quad_data_norm;
lap_vec_custom = zeros(size(quad_data_norm));

for i = 1 : (length(IR_ind) - 1)
    curr_ind = IR_ind(i) - 2 : IR_ind(i + 1) - 3;
%     I have times when the rew is with the IR. no double inds
%     allowed this way, plus jitter around IR allowed. 
    
    if sum(ismember(curr_ind, valve_ind))~= 1
       disp(['bad lap ', num2str(i)])
       remove_us_vec(curr_ind) = 1;
       quad_data(find(remove_us_vec)) = nan;
    end
    
    curr_quad = quad_data(curr_ind);
    
    if sum(isnan(curr_quad)) > 0 && ~isnan(curr_quad(1))
        curr_quad(1) = nan;
        min_quad_ind = find(~isnan(curr_quad), 1);
    else
        min_quad_ind = 1;
    end
    quad_data_norm(curr_ind) = ...
        mod((curr_quad - curr_quad(min_quad_ind))./...
        (curr_quad(end) - curr_quad(min_quad_ind)) + IR_shift, 1); 
    
    % such that the IR is in desired bin for comparing to other trials
    
    discrete_quad_data(curr_ind) = ...
        discretize(mod(quad_data_norm(curr_ind), max(quad_data_norm(curr_ind))), ...
        linspace(0, 1, bin_number + 1));
    lap_vec_custom(IR_ind(i) - 2) = 1;
end
lap_vec_custom(IR_ind(end)) = 1;

lap_vec_custom = cumsum(lap_vec_custom, 'omitnan');

% lap division - according to the IR, not the rewards!!

%% visualize lap devision:
figure;
subplot(1,2,1)
scatter(1 : length(quad_data_norm), quad_data_norm, 5,...
    discrete_quad_data, 'filled');
colormap(gca, 'hsv')
hold on;
plot(IR_ind, quad_data_norm(IR_ind), '^k');
plot(valve_ind, quad_data_norm(valve_ind), 'kX');
legend('path','IR', 'valve')

subplot(1,2,2)
scatter(1 : length(quad_data_norm), quad_data_norm, 5,...
    lap_vec_custom, 'filled');
colormap(gca, 'prism')

end