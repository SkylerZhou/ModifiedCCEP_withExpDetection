function plot_ccep_pair(ax, out, stim_ch, resp_ch, line_color, n1_time, n2_time, zoom_times, zoom_factor, which, nonsig_overlap_idx)

%% extract data corresponding to the specific stim and resp
row = find(strcmp(out.bipolar_labels, stim_ch));
col = find(strcmp(out.bipolar_labels, resp_ch));
which_n = 1;
    
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




%% get why it was rejected for the nonsig signals 
if ~isnan(nonsig_overlap_idx) % if we have started to plot the nonsig signals 
    sig_avg = out.rejection_details(which_n).reject.sig_avg;
    pre_thresh = out.rejection_details(which_n).reject.pre_thresh;
    at_thresh = out.rejection_details(which_n).reject.at_thresh;
    %no_both = out.rejection_details(which_n).reject.no_both;
    exp = out.rejection_details(which_n).reject.exp;

    why = nan;

    if sig_avg(row,col) == 1
        why = 'averaging';
    end
    if pre_thresh(row,col) == 1
        why = 'artifact';
    end
    if at_thresh(row,col) == 1      
        why = 'threshold';
    end
    if exp(row,col) == 1
        if isnan(why)
            why = 'exponential';
        end
    end
end
    %{
    %if no_both(row,col) == 1
    %    if isnan(why)
    %        why = 'no both';
    %    end
    %end
    %}


%% for text 
stim_label = num2str(stim_ch);
resp_label = num2str(resp_ch);
stim_start = out.other.periods(row).start_time;

%% plot 
plot(ax, eeg_times, avg, 'Color', line_color, 'LineWidth', 1.5);
hold(ax, 'on');

% N1 annotation
if ~isnan(out.elecs(row).N1(col, 1)) 
    x = out.elecs(row).N1(col, 2) / out.other.stim.fs;
    x_indx = round(out.elecs(row).N1(col, 2) + stim_idx + 1);
    y = out.elecs(row).detrend_filt_avgs(x_indx, col);
    plot(ax, x, y, 'cX', 'MarkerSize', 4, 'LineWidth', 1.5);
    text(ax, wav_time + 0.01, avg(round(wav_idx)), sprintf('%s z-score: %1.1f', which, wav(1)), 'FontSize', 6);
end

% N2 annotation
if ~isnan(out.elecs(row).N2(col, 2))
    x = out.elecs(row).N2(col, 2) / out.other.stim.fs;
    x_indx = out.elecs(row).N2(col, 2) + stim_idx + 1;
    y = out.elecs(row).detrend_filt_avgs(x_indx, col);
    plot(ax, x, y, 'mX', 'MarkerSize', 4, 'LineWidth', 1.5);
end

% set xlim and ylim 
xlim(ax, [zoom_times(1) zoom_times(2)]);
height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2)) - median(avg)));
if ~any(isnan(avg))
    ylim(ax, [median(avg) - zoom_factor * height, median(avg) + zoom_factor * height]);
end

% vline for time=0
plot(ax, [0 0], ylim(ax), 'k--');

% for the nonsig subplots, specify reasons of rejection
if ~isnan(nonsig_overlap_idx)
    xl = xlim;
    yl = ylim;
    text(ax, xl(1),yl(2), sprintf('Stim: %s\nResponse: %s\nStart: %.2f s\nRejected: %s',stim_label,resp_label,stim_start,why),...
        'FontSize', 5, ...
        'horizontalalignment','left',...
        'verticalalignment','top');
else
    xl = xlim;
    yl = ylim;
    text(ax, xl(1),yl(2), sprintf('Stim: %s\nResponse: %s\nStart: %.2f s',stim_label,resp_label,stim_start),...
        'FontSize', 5, ...
        'horizontalalignment','left',...
        'verticalalignment','top');
end

end