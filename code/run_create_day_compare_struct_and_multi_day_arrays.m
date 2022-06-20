
%% create structure with all the files for all the days
data_folder= 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_folder);
file_names = dir('*_full_sig_ROIs_df_F.signals');
file_names = {file_names.name};
mouse_names = ...
    cellfun(@(x) (x(1:3)), file_names, 'uniformoutput', false);
[~, mouse_names] = findgroups(mouse_names);

for i = 1 : length(mouse_names)
        
    day_compare_struct = struct();
    
    if strcmp(mouse_names{i}, '8C3') || strcmp(mouse_names{i}, '9B2')
        day_code = {'100', '000', '001', '002'}; 
    else
        day_code = ...
            {'000', '001', '002', '003', '004', '005'}; 
    end
    
    for j = 1 : length(day_code)
        day_compare_struct.(['day_' num2str(day_code{j})]) = [];
    end
    
    curr_file_names = dir([mouse_names{i} '*_full_sig_ROIs_df_F.signals']);
    curr_file_names = {curr_file_names.name}';
    ind = cellfun(@(x) strfind(x,'_'), curr_file_names, 'uniformoutput', false);
    ind = cellfun(@(x) x(3)-1, ind, 'uniformoutput', false);
    
    curr_file_names = ....
        cellfun(@(x, y) x(1 : y), curr_file_names, ind, 'uniformoutput', false);

    
    day_compare_struct = ...
        create_day_compare_struct_NEW(day_compare_struct, curr_file_names, day_code);
    
    save([mouse_names{i}, '_day_compare_struct.mat'], 'day_compare_struct');
    
    %% create and save all the required variables for comparing days
    % the variables saved here: sig_array, event_array, mask_struct
    day_fields = fieldnames(day_compare_struct);
    
    mask_struct = create_multi_day_arrays(day_compare_struct, day_fields);
    save([mouse_names{i}, '_mask_struct_df_F_all_ROIs.mat'], 'mask_struct');
    
    
end

function day_compare_struct = ...
    create_day_compare_struct_NEW(day_compare_struct, file_names, day_code)

day_fields = fieldnames(day_compare_struct);
%

for i = 1 : length(file_names)
    curr_session = file_names{i};
    
    curr_day = strfind(day_code, curr_session(5:7));
    
    curr_day = find(~cellfun(@isempty, curr_day));
    
    day_compare_struct.(day_fields{curr_day}) = curr_session;
end
end