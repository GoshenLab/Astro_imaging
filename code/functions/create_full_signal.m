function [full_sig, ot_length, cond_full] = create_full_signal(file_name,...
    frame_num, start_here, corr_ind, varargin)

% corr_ind is the red channel based movement detection for 2 ot's (if both
% are true, should be included).
% varagin - if there's additional elimination according to
% elimination_movie_based...

disp('creating full sig');

%% create full signal - omitting the first frames:
sig_file = dir([file_name, '_ot_000*.signals']);
if ~isempty(sig_file)
    load(sig_file.name, '-mat')
    if ndims(sig) == 3
        sig = squeeze(sig(1, :, :));
    end
        sig_interp.ot_000 = ...
            interp1(1 : 2 : frame_num, sig, 1 : frame_num);  
else 
       sig_interp.ot_000 = [];
end

sig_file = dir([file_name, '_ot_001*.signals']);
if ~isempty(sig_file)
    load(sig_file.name, '-mat')
    if ndims(sig) == 3
        sig = squeeze(sig(1, :, :));
    end
else
    sig = [];
end
sig_interp.ot_001 = interp1(2 : 2 : frame_num, sig, 1 : frame_num);

full_sig = [sig_interp.ot_000(start_here : end, :),...
    sig_interp.ot_001(start_here : end, :)];
ot_length = [size(sig_interp.ot_000, 2), size(sig_interp.ot_001, 2)];
cond_full = sum(corr_ind, 2) == 2;
save([file_name, '_full_sig_V2.mat'], 'full_sig', 'ot_length');

if ~isempty(varargin)
    bad_frames_from_movie = varargin{1};
    cond_full(bad_frames_from_movie) = false;
end


end
