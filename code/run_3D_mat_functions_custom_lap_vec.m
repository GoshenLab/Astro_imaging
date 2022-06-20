data_folder_name = 'D:\astro_imaging\Nature_code\all_code_data\';
cd(data_folder_name);

mouse_names = {...
    '9Q4',...
    '8C3',...
    '9B2',...
    '9H3',...
    '9I3',...
    '9P2',...
    '9N1',...
    '9O3',...
    '9S3',...
    };

%%
for m = 1 : length(mouse_names)
    mouse_name = mouse_names{m};
    load([mouse_name, '_mask_struct_df_F_all_ROIs.mat'])
    
    file_names = {mask_struct.name};
    how_many_cells = [mask_struct.how_many_cells];
    
    for discrete_analysis = [false] % can run on analouge too if needed
        for j = 1 : length(file_names)
            
            file_name = file_names{j};
            
            load([file_name, '_full_sig_ROIs_df_F.signals'], '-mat')
            load([file_name, '_new_event_det_df_F.mat'], 'events_above_min');
            event_mat = double(events_above_min);
            
            load([file_name, '_remove_us_vec_trimmed.mat'])
            remove_us_vec = logical(remove_us_vec);
            df_F(remove_us_vec, :) = nan;
            event_mat(remove_us_vec, :) = nan;
            
            if discrete_analysis
                save_name = [file_name, '_discrete'];
                sig = event_mat;
            else
                save_name = [file_name, '_analoge'];
                sig = df_F;
            end
            
            disp(save_name);
            
            all_ROIs = 1 : how_many_cells(j);
            
            load([file_name,...
                '_lap_vec_custom_discrete_quad_data_quad_data_norm.mat'],...
                'lap_vec_custom',...
                'discrete_quad_data');
            lap_vec = lap_vec_custom;
            load([file_name,...
                '_quad_data_shift_movement_velocity_hardware.mat'],...
                'movement');
            
            movement = logical(movement);
           
            good_frames = ...
                (sum(isnan(sig), 2) == 0) & (movement == 1);
            
            %% location 3D mat
            if ~isfile([save_name, '_sig_3D_mat_df_F_custom_lap.mat'])
                disp('calculating loc_3D_mat')
                sig_3D_mat = ...
                    create_3d_activity_mat(lap_vec(good_frames),...
                    discrete_quad_data(good_frames),...
                    sig(good_frames, :));
                disp('saving loc_3D_mat')
                save([save_name, '_sig_3D_mat_df_F_custom_lap'], 'sig_3D_mat');
            else
               disp('loc_3D_mat already created') 
            end

        end
    end
end
