function title_str = plot_ccep_pair(ax, out, stim_ch, resp_ch, line_color, n1_time, n2_time, zoom_times, zoom_factor, which)

% extract data corresponding to the specific stim and resp
row = find(strcmp(out.bipolar_labels, stim_ch));
col = find(strcmp(out.bipolar_labels, resp_ch));
    
avg = out.elecs(row).detrend_filt_avgs(:, col);
times = out.elecs(row).times;
eeg_times = convert_indices_to_times(1:length(avg), out.other.stim.fs, times(1));
wav = out.elecs(row).(which)(col, :);
                
stim_idx = out.elecs(row).stim_idx;
wav_idx = wav(2) + stim_idx + 1;
wav_time = convert_indices_to_times(wav_idx, out.other.stim.fs, times(1));
n1_idx = floor(n1_time * out.other.stim.fs);
n2_idx = floor(n2_time * out.other.stim.fs);
temp_n1_idx = n1_idx + stim_idx - 1;

% plot 
plot(ax, eeg_times, avg, 'Color', line_color, 'LineWidth', 1.5);
title_str = sprintf('Stim %s, Resp %s', num2str(stim_ch), num2str(resp_ch));
hold(ax, 'on');

% N1 annotation
if (out.elecs(row).N1(col, 1)) ~= 0
    x = out.elecs(row).N1(col, 2) / out.other.stim.fs;
    x_indx = round(out.elecs(row).N1(col, 2) + stim_idx + 1);
    y = out.elecs(row).detrend_filt_avgs(x_indx, col);
    plot(ax, x, y, 'bX', 'MarkerSize', 8, 'LineWidth', 1.5);
    text(ax, wav_time + 0.01, avg(round(wav_idx)), sprintf('%s z-score: %1.1f', which, wav(1)), 'FontSize', 6);
end

% N2 annotation
if ~isnan(out.elecs(row).N2(col, 2))
    x = out.elecs(row).N2(col, 2) / out.other.stim.fs;
    x_indx = out.elecs(row).N2(col, 2) + stim_idx + 1;
    y = out.elecs(row).detrend_filt_avgs(x_indx, col);
    plot(ax, x, y, 'rX', 'MarkerSize', 8, 'LineWidth', 1.5);
end

% set xlim and ylim 
xlim(ax, [zoom_times(1) zoom_times(2)]);
height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2)) - median(avg)));
if ~any(isnan(avg))
    ylim(ax, [median(avg) - zoom_factor * height, median(avg) + zoom_factor * height]);
end

% vline for time=0
plot(ax, [0 0], ylim(ax), 'k--');

end