function [df_F,  deltaF_dombeck_sig_baseline] =...
    calculate_df_F( sig, time_win, sampling_rate, per_centile)
% sig - can be matrix
% time_win - take X sec window around point
% percentile - to smooth according to

sig_length = size(sig, 1);
sig_number = size(sig, 2);
sample_win = round(time_win * sampling_rate); % sample_win - number of samples

if mod(sample_win,2) == 1
    sample_win = sample_win - 1;
end


%% correct slow changes in F

sig_pad = ...
    [nan(0.5 * sample_win, sig_number);...
    sig; ...
    nan(0.5 * sample_win, sig_number)];

% create matrix for 'convolution'
x = [1 : sig_length]';
x = repmat(x, 1, sample_win);
x2 = [1 : sample_win] - 1;
x2 = repmat(x2, sig_length, 1);
X = x + x2;
deltaF_dombeck_sig_baseline = zeros(size(sig));

for s = 1 : size(sig, 2)
    disp(s)
    temp = sig_pad(:, s);
    deltaF_dombeck_sig_baseline(:, s) = prctile((temp(X))', per_centile);
end

%%
df_F = (sig - deltaF_dombeck_sig_baseline) ./ deltaF_dombeck_sig_baseline;
end