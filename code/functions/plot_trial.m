function h_f = plot_trial(frame_num, frame_rate, belt_length, back_trace,...
    valve_start, lickometer_ind, lap_ind_IR, remove_us_vec)

% back trace - the trace on which the hardware activity would be indicated

h_f = figure;
colormap hsv;
p0 = scatter((1 : frame_num) ./ (frame_rate * 60), back_trace, 15,...
    back_trace./max(back_trace), 'filled');
hold on
p2 = plot(lickometer_ind / (frame_rate * 60), back_trace(lickometer_ind),...
    'o', 'MarkerSize', 3, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
pIR = plot(lap_ind_IR / (frame_rate * 60), back_trace(lap_ind_IR),...
    'k^', 'MarkerSize', 10, 'LineWidth', 2);
p1 = plot(valve_start / (frame_rate * 60), back_trace(valve_start),'kX',...
    'MarkerSize', 20, 'LineWidth', 2);
xlabel('time (min)', 'FontSize', 18)
ylabel('distance (cm)','FontSize', 18)
xlim([0 ceil(frame_num / (frame_rate * 60))])
ylim([0 belt_length])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

plot(find(remove_us_vec) ./ (frame_rate * 60), back_trace(find(remove_us_vec)),...
    'y.', 'MarkerSize', 15);

legend('path', 'lick', 'IR', 'valve')
end
