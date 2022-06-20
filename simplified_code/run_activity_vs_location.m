
ALL_DATA = struct();
ALL_DATA_matching = struct();
only_matching = true;   % test all the ROIs or only the repeated active ones
lap_cutoff = 20; % can be [] if you want all laps.
quad_or_vel = 1;


SAMPLING_RATE = 15.49;

BELT_LENGTH = 170;
PERM_NUM = 1000;

bin_number_vel = 10;
bin_number_quad = 10;

discrete_analysis = true; % false for df/f

% include_cells - active cells (and repeated if needed)
                  
              
[template_mat, template_mat_perm] = ...
    deal_with_file(include_cells, lap_cutoff,...
    quad_or_vel, discrete_analysis, PERM_NUM,...
    bin_number_vel, bin_number_quad);

mean_for_each_perm = ...
    squeeze(mean(template_mat_perm, 2, 'omitnan'));
mean_across_all_perm_V1 = mean(mean_for_each_perm, 2, 'omitnan');
            

plot_template_mat_results(template_mat, mean_across_all_perm_V1)
            
if  quad_or_vel == 1
    x = 2 : 10;         % omit the reward bin (bin 1)
else
    x = 1 : 10;         % this is also the case when using random rewards
end

template_mat = template_mat(x, :);
filler = repmat(x, size(template_mat, 2), 1)';

R_real = corr(filler(:), template_mat(:), 'rows', 'pairwise');
R_perm = nan(size(template_mat_perm, 3), 1);
for s = 1 : size(template_mat_perm, 3)
    curr_data = template_mat_perm(x, :, s);
    curr_data = curr_data(:);
    R_perm(s) = corr(filler(:), curr_data(:), 'rows', 'pairwise');
end

R_real_p_val = ...
    sum(R_perm > R_real)./PERM_NUM;



%% functions:

function [template_mat, template_mat_perm] = ...
    deal_with_file(lap_cutoff, quad_or_vel,...
    discrete_analysis, PERM_NUM, bin_number_vel, bin_number_quad)

% load the relevant files (lap_vec, df_F, event_mat, discrete_location,...
% movement, raw_velocity)

if ~isnan(lap_cutoff)
    cut_point = find(lap_vec == (lap_cutoff + 1), 1);
    if ~isempty(cut_point)
        df_F(cut_point:end, :) = nan;
        event_mat(cut_point:end, :) = nan;
        % to avoid non-existing values
        velocity(cut_point:end, :) = nan;
    end
end

movement = logical(movement);

velocity(~movement) = nan;
velocity(...
    isnan(discrete_location)) = nan; 
vel_data_norm = ...
    normalize(velocity, 'range');
discrete_vel_data = ...
    discretize(vel_data_norm, bin_number_vel);

%% I left some constants here if I need to plot this.
if quad_or_vel == 1
    discrete_data = discrete_location;
    bin_num = bin_number_quad;
else
    discrete_data = discrete_vel_data;
    bin_num = bin_number_vel;
end

if discrete_analysis
    sig = event_mat;
else
    sig = df_F;    
end
%%

if discrete_analysis
    concurrent_event_num = sum(sig, 2);
else
    concurrent_event_num = mean(sig, 2);
end


good_frame_ind = ...
    ~((sum(isnan(sig), 2) > 0) |...
    (movement ~= 1));


template_mat = ...
    number_of_sim_events_in_bin...
    (lap_vec(good_frame_ind),...
    discrete_data(good_frame_ind),...
    concurrent_event_num(good_frame_ind), 1);


if discrete_analysis
    % this is in order to get the cell proportion!
    template_mat = template_mat./size(event_mat, 2);
end

%% permutation
sig_for_perm = sig(good_frame_ind, :);

perm_delta = randperm(length(sig_for_perm), min(PERM_NUM, size(sig_for_perm, 1)));
disp('calculating template_mat_perm')
template_mat_perm = ...
    nan(bin_num, size(template_mat, 2), min(PERM_NUM, size(sig_for_perm, 1)));

for p = 1 : min(PERM_NUM, size(sig_for_perm, 1))
    
    concurrent_event_num_perm = ...
        circshift(concurrent_event_num(good_frame_ind),...
        perm_delta(p));
    
    temp = number_of_sim_events_in_bin...
        (lap_vec(good_frame_ind),...
        discrete_data(good_frame_ind),...
        concurrent_event_num_perm, 1);
    %
    %% normalization for discrete analysis:
    if discrete_analysis
        temp = temp./size(event_mat, 2);
    end
    %
    template_mat_perm(1 : size(temp, 1),...
        1 : size(temp, 2), p) = ...
        temp;
end



end

function plot_template_mat_results(template_mat, mean_across_all_perm_V1)

x = [1 : size(template_mat, 1)];
shadedErrorBar(x,...
    mean(template_mat, 2, 'omitnan') ./ mean_across_all_perm_V1,...
    std(template_mat, [], 2, 'omitnan')./...
    sqrt(sum(~isnan(template_mat),2))./...
    mean_across_all_perm_V1);
hold on;
plot([min(x) max(x)], [1 1], 'k--')

end