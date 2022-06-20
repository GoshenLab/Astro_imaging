

data_dir = 'D:\astro_imaging\Nature_code\all_code_data\random_reward_data\';
cd(data_dir)

%% create df/f without noisy frames:
time_win = 250;          % slow astro CA2 signal
sampling_rate = 15.49/2; % using the optotune
per_centile = 8;

deltaF_dombeckMat_param = struct();
deltaF_dombeckMat_param.time_win = time_win;
deltaF_dombeckMat_param.sampling_rate = sampling_rate;
deltaF_dombeckMat_param.per_centile = per_centile;

all_file_names = dir('*_quadrature.mat');
all_quad_file_name = {all_file_names.name};
file_names = cellfun(@(x) x(1:11), all_quad_file_name, 'uniformoutput', false);



for i = 1 : length(file_names)
    disp( file_names{i} );
    if isempty(dir([file_names{i}, '_ot_001*.signals'])) & ...
            isempty(dir([file_names{i}, '_ot_000*.signals']))
        disp('no signals on both ots')
        continue
    end
    
    if ~isfile([file_names{i} '_full_sig_ROIs_df_F.signals'])
        
    load([file_names{i} '_frame_num.mat'])
    pipeline_revised_func(file_names{i}, frame_num,...
        deltaF_dombeckMat_param);
    else
        disp(['already created df/F for file ', file_names{i}]);
    end
end

%%


function pipeline_revised_func(file_name, frame_num, deltaF_dombeckMat_param)

start_here = 6;

if isfile([file_name, '_corr_ind.mat'])
    load([file_name, '_corr_ind.mat'], 'corr_ind')
else
    corr_ind = ones(frame_num - start_here + 1, 2);
end


if isfile([file_name, '_elimination_movie_based.mat'])
    load([file_name, '_elimination_movie_based.mat'],...
        'elimination_movie_based');
else
    elimination_movie_based.start = [];
    elimination_movie_based.end = [];
end

ind_vec = 1 : frame_num;
if isempty(elimination_movie_based.start) & ...
        isempty(elimination_movie_based.end)
    elimination_movie_based.start = 1;
    elimination_movie_based.end = frame_num;
else
    disp('elimination by movie is relevant')
    if isempty(elimination_movie_based.start)
        elimination_movie_based.start = 1;
    end
    if isempty(elimination_movie_based.end)
        elimination_movie_based.end = frame_num;
    end
    
end

bad_frames_from_movie = ...
    (ind_vec < elimination_movie_based.start |...
    ind_vec > elimination_movie_based.end);
bad_frames_from_movie = bad_frames_from_movie(start_here : end);

if ~isfile([file_name, '_full_sig_V2.mat'])
    [full_sig, ot_length, cond_full] = create_full_signal( file_name,...
        frame_num, start_here, corr_ind, bad_frames_from_movie); 
else
    load([file_name, '_full_sig_V2.mat'])
    cond_full = sum(corr_ind, 2) == 2;
    cond_full(bad_frames_from_movie) = false;
end

%% this is for visualizing omitted frames
figure;
if ~isempty(bad_frames_from_movie)
    cond_full(bad_frames_from_movie) = false;
end
x = 1 : length(cond_full);
plot(x, full_sig + [1:size(full_sig, 2)]*10^4, 'k');
hold on;
plot(x(cond_full), full_sig(cond_full,:) + [1:size(full_sig, 2)]*10^4, '.');
%%

full_sig_for_df = full_sig;
full_sig_for_df(~cond_full, :) = nan;

deltaF_dombeckMat_param.ot_length = ot_length;

df_F =...
    calculate_df_F( full_sig_for_df, deltaF_dombeckMat_param.time_win,...
    deltaF_dombeckMat_param.sampling_rate,...
    deltaF_dombeckMat_param.per_centile);


save([file_name '_full_sig_ROIs_df_F.signals'],...
    'df_F', 'deltaF_dombeckMat_param');

end



