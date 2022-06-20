function index = find_events_with_large_ISI(data, threshold_sig, min_ISI)

a = data > threshold_sig;
[~,index] = findpeaks(double(a), 'MinPeakDistance', min_ISI);

end