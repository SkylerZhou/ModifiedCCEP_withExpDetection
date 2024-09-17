% load ori_out and new_out
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_idx = 45;
which_n = 1;

ori_patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'ori_pipeline', patient_files(patient_idx));
temp = load(ori_patient_file);
ori_out = temp.pt_out;
    
new_patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(patient_idx));
temp = load(new_patient_file);
new_out = temp.new_out;

% get rejection info
sig_avg_ori = ori_out.rejection_details(which_n).reject.sig_avg;
pre_thresh_ori = ori_out.rejection_details(which_n).reject.pre_thresh;
at_thresh_ori = ori_out.rejection_details(which_n).reject.at_thresh;

any_reject_ori = sig_avg_ori == 1 | pre_thresh_ori == 1 | at_thresh_ori == 1;
exp_new = new_out.rejection_details(which_n).reject.exp;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For those rejected by 'exp', randomly select and plot 25 figures
% to see how good curve fitting gof is in rejecting true artifacts 

% Parameters
number_to_plot = 25;
zoom_times = [-300e-3 300e-3];
peak_end_time = 0.3;
aLower = -5;
aUpper = 5;
bLower = -50;
bUpper = 0;

% Find row and col for waves that were only rejected due to 'exp' in
% 'new_out' but not rejected in 'original_out'
exp_reject = (exp_new == 1); 
%[row, col] = find(exp_reject);
[row, col] = find(exp_reject & ~any_reject_ori);

% Randomly select pairs
rand_indices = randperm(length(row), number_to_plot);

% Initialize first figure for after_peak_avg_norm plots
figure(1);
set(gcf, 'position', [100, 100, 1200, 1000]);
t1 = tiledlayout(5, 5, 'padding', 'compact', 'tilespacing', 'compact');
title(t1, 'Standarized After-Peak Avg with Exponential Fit');

% Initialize second figure for avg plots
figure(2);
set(gcf, 'position', [100, 100, 1200, 1000]);
t2 = tiledlayout(5, 5, 'padding', 'compact', 'tilespacing', 'compact');
title(t2, 'Original Avg with Exponential Fit');

% Loop to plot random selections
for i = 1:number_to_plot
    stim = row(rand_indices(i));
    resp = col(rand_indices(i));

    % Get labels 
    stim_label = new_out.bipolar_labels{stim};
    resp_label = new_out.bipolar_labels{resp};

    % Get the average waveform (avg) and corresponding time values (eeg_times)
    avg = new_out.elecs(stim).detrend_filt_avgs(:, resp);
    times = new_out.elecs(stim).times;
    eeg_times = convert_indices_to_times(1:length(avg), new_out.other.stim.fs, times(1));

    % Plot gof and fit the exponential function
    stim_idx = new_out.elecs(stim).stim_idx;
    peak_start_index = new_out.elecs(stim).N1(resp,2) + stim_idx;  % the index of 0 sec is stim_idx
    peak_end_index = stim_idx + floor(new_out.other.stim.fs * peak_end_time);                  
    peak_times = convert_indices_to_times(peak_start_index:peak_end_index, new_out.other.stim.fs, times(1));
    after_peak_avg = new_out.elecs(stim).detrend_filt_avgs(peak_start_index:peak_end_index, resp);                 
    
    % standarize the after_peak_avg 
    mean_val = mean(after_peak_avg);
    std_val = std(after_peak_avg);
    if std_val ~= 0
        after_peak_avg_std = (after_peak_avg - mean_val) / std_val;
    else
        after_peak_avg_std = zeros(size(after_peak_avg)); 
    end

    % Define cutoff frequency and filter parameters
    cutoff_freq = 50;  
    Fs = 2 * new_out.other.stim.fs;
    [b, a] = butter(4, cutoff_freq/(Fs/2), 'low');  % Fs is the sampling frequency
    after_peak_avg_std = filtfilt(b, a, after_peak_avg_std);

    % Examine goodness of fit using exponential fit
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
                            'Lower', [aLower, bLower], ...
                            'Upper', [aUpper, bUpper]);
    [f, gof] = fit(peak_times', after_peak_avg_std(:), 'exp1', fitOptions);
    a = f.a;
    b = f.b;

    % Diagnostic checks
    fprintf('Iteration %d: stim = %d, resp = %d\n', i, stim, resp);
    fprintf('Size of avg: %s\n', mat2str(size(avg)));
    fprintf('Size of eeg_times: %s\n', mat2str(size(eeg_times)));

    if isempty(avg)
        warning('avg is empty for stim %d and resp %d.', stim, resp);
        continue;
    end

    if all(isnan(avg))
        warning('avg contains only NaNs for stim %d and resp %d.', stim, resp);
        continue;
    end

    if length(avg) ~= length(eeg_times)
        warning('Mismatch in lengths of avg and eeg_times for stim %d and resp %d.', stim, resp);
        continue;
    end

    %% Plot after_peak_avg_norm and the fitted exponential function in Figure 1
    ax1 = nexttile(t1);
    plot(ax1, peak_times, after_peak_avg_std, 'k', 'linewidth', 2); % Plot std waveform
    hold(ax1, 'on');

    % Fitted exponential function
    plot(ax1, peak_times, f(peak_times), 'r--', 'linewidth', 1.5);  % Fitted curve in red dashed line

    % Set the axis limits
    xlim(ax1, zoom_times);
    ylim(ax1, [-1.2 * max(abs(after_peak_avg_std)), 1.2 * max(abs(after_peak_avg_std))]);

    % Place the annotation
    xl = xlim(ax1);
    yl = ylim(ax1);
    text(ax1, xl(1) + 0.02 * diff(xl), yl(2) - 0.05 * diff(yl), ...
        sprintf('Stim: %s\nResp: %s\nGOF: %.2f\na: %.2f\nb: %.2f', stim_label, resp_label, gof.adjrsquare, a, b), ...
        'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', 9);

    % Mark the stimulation point
    plot(ax1, [0 0], ylim(ax1), 'k--');

    %% Plot avg and the fitted exponential function in Figure 2
    ax2 = nexttile(t2);
    plot(ax2, eeg_times, avg, 'k', 'linewidth', 2); % Plot original average waveform
    hold(ax2, 'on');

    % Fitted exponential function
    plot(ax2, peak_times, f(peak_times), 'r--', 'linewidth', 1.5);  % Fitted curve in red dashed line

    % Set the axis limits
    xlim(ax2, zoom_times);
    ylim(ax2, [-1.2 * max(abs(avg)), 1.2 * max(abs(avg))]);

    % Place the annotation
    xl = xlim(ax2);
    yl = ylim(ax2);
    text(ax2, xl(1) + 0.02 * diff(xl), yl(2) - 0.05 * diff(yl), ...
        sprintf('Stim: %s\nResp: %s\nGOF: %.2f\na: %.2f\nb: %.2f', stim_label, resp_label, gof.adjrsquare, a, b), ...
        'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', 9);

    % Mark the stimulation point
    plot(ax2, [0 0], ylim(ax2), 'k--');
end
