function [sig_3D_mat, G, template_mat] = ...
    create_3d_activity_mat(lap_vec, discrete_data, sig)
% INPUT:
%   lap_vec:        vector indicating lap vector (horizontal).
%                   if not nessacry, should be ones. the lap_vec_custom
%                   starts at 0, but since this lap is not complete, the
%                   quad_data_discrete will include NANs at the begining.
%                   Careful if you use a different variable.
%   discrete_data:  discretized data (e.g. quad_data), horizontal.
%                  
%   sig:            the relevant signal

% OUTPUT
%   sig_3D_mat:     discrete_data bins X lap X ROI

% find the groups in the data, create table
group_table = table();

group_table.lap_vec = ...
    reshape(lap_vec, length(lap_vec), 1);
group_table.discrete_data = ...
    reshape(discrete_data, length(discrete_data), 1);
G = findgroups(group_table);

% don't include NAN data, prepare the data template per ROI
bin_number = max(group_table.discrete_data);
include_ind = ~isnan(G);
template_mat = nan(max(group_table.discrete_data(include_ind)),...
    max(lap_vec(include_ind)));
template_ind = group_table.discrete_data + (group_table.lap_vec - 1) * bin_number;
template_ind = template_ind(~isnan(template_ind));
template_mat(template_ind) = 1;

% intiate the 3d mat:
sig_3D_mat = zeros(max(group_table.discrete_data(include_ind)),...
    max(group_table.lap_vec(include_ind)),...
    size(sig, 2));

% insert the ROI data into the 3d mat:
for i = 1 : size(sig, 2)
    curr_data = splitapply(@(x) mean(x, 'omitnan'),...
        sig(include_ind, i), G(include_ind));
    temp = template_mat;
    temp(temp == 1) = curr_data;
    sig_3D_mat(:, :, i) = temp;
    
end
end