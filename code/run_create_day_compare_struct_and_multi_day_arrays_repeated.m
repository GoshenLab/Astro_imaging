data_folder= 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_folder);

mouse_name = {...
'9Q4';    
'8C3';
    '9B2';
    '9P2';
    '9N1';
    '9N3';
    '9O3';
    '9O2';
    '9S3';
    };

repeated_days = [...
3;    
2;
    2;
    3;
    3;
    2;
    2;
    2;
    3];

for i = 1 : length(mouse_name)
   load([mouse_name{i} '_mask_struct_df_F_all_ROIs.mat'])
   mask_struct = mask_struct(1 : repeated_days(i));
   save([mouse_name{i} '_mask_struct_df_F_matching_ROIs.mat'],...
       'mask_struct')
   
end