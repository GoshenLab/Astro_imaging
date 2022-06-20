function template_mat = number_of_sim_events_in_bin(discrete_var1,...
    discrete_var2, sig, mean_or_sum)

% create matrix of mean signal (mean num of concurrent events/df_f) using 2
% discrete grouping variables (e.g. lap number, location, velocity)
% sig is a vector (population mean/single ROI) 

discrete_var1 = reshape(discrete_var1, length(discrete_var1), 1);
discrete_var2 = reshape(discrete_var2, length(discrete_var2), 1);

if sum(discrete_var1 == 0) > 0
    discrete_var1 = discrete_var1 + 1;
end
if sum(discrete_var2 == 0) > 0
    discrete_var2 = discrete_var2 + 1;
end

group_table = table();
group_table.discrete_var1 = discrete_var1;
group_table.discrete_var2 = discrete_var2;
G = findgroups(group_table);


% don't include NAN data, prepare the data template
include_ind = ~isnan(G);
template_mat = nan(max(discrete_var2(include_ind)),...
    max(discrete_var1(include_ind)));
template_ind = discrete_var2 + (discrete_var1-1) * max(discrete_var2);
template_ind = template_ind(~isnan(template_ind));
template_mat(template_ind) = 1;

if mean_or_sum == 1
    sig = reshape(sig, length(sig), 1);
    temp = splitapply(@(x) mean(x, 'omitnan'),...
        sig, G);    
elseif mean_or_sum == 2
    % this is to see how many bins are in this condition
    temp = splitapply(@(x) sum(x, 'omitnan'),...
        sig, G);
end
template_mat(template_mat == 1) = temp;

end