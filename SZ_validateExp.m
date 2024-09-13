% load ori_out and new_out
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_idx = 40;
which_n = 1;

ori_patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(patient_idx));
temp = load(ori_patient_file);
ori_out = temp.pt_out;
    
new_patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'modified_pipeline', patient_files(patient_idx));
temp = load(new_patient_file);
new_out = temp.new_out;

% get rejection info
sig_avg_ori = ori_out.rejection_details(which_n).reject.sig_avg;
pre_thresh_ori = ori_out.rejection_details(which_n).reject.pre_thresh;
at_thresh_ori = ori_out.rejection_details(which_n).reject.at_thresh;

any_reject_ori = sig_avg_ori == 1 | pre_thresh_ori == 1 | at_thresh_ori == 1;
exp_new = new_out.rejection_details(which_n).reject.exp;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for those rejected by exp, randomly select and plot 25 figures
%% to see how good curve fitting gof is in rejecting true artifacts 

% Parameters
number_to_plot = 25;
zoom_times = [-300e-3 300e-3];
peak_end_time = 0.3;
aLower = -2000;
aUpper = 2000;
bLower = -50;
bUpper = 0;


% find row and col for waves that were only rejected due to 'exp' in
% 'new_out' but not rejected in 'original_out'
exp_reject = (exp_new == 1); 
[row,col] = find(exp_reject & ~any_reject_ori);


% Randomly select pairs
rand_indices = randperm(length(row), number_to_plot);


% Initialize figure
figure;
set(gcf, 'position', [100, 100, 1200, 1000]);
t = tiledlayout(5, 5, 'padding', 'compact', 'tilespacing', 'compact');


% Loop to plot random selections
for i = 1:number_to_plot
    stim = row(rand_indices(i));
    resp = col(rand_indices(i));

    % Get the average waveform (avg) and corresponding time values (eeg_times)
    avg = new_out.elecs(stim).detrend_filt_avgs(:, resp);
    times = new_out.elecs(stim).times;
    eeg_times = convert_indices_to_times(1:length(avg), new_out.other.stim.fs, times(1));

    % Plot gof and fit the exponential function
    stim_idx = new_out.elecs(stim).stim_idx;
    peak_start_index = new_out.elecs(stim).N1(resp,2) + stim_idx;  % the index of 0 sec is 512
    peak_end_index = stim_idx + floor(new_out.other.stim.fs * peak_end_time);                  
    peak_times = convert_indices_to_times(peak_start_index:peak_end_index, new_out.other.stim.fs, times(1));
    after_peak_avg = new_out.elecs(stim).detrend_filt_avgs(peak_start_index:peak_end_index, resp);                 
    
    % Examine goodness of fit using exponential fit
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'Lower', [aLower, bLower], ...
                        'Upper', [aUpper, bUpper]);
    [f, gof] = fit(peak_times', after_peak_avg(:), 'exp1', fitOptions);
    a = f.a;
    b = f.b;

    % Plot in the next tile
    nexttile;
    plot(eeg_times, avg, 'k', 'linewidth', 2); % Plot original average waveform
    hold on;

    % Plot the fitted exponential function
    plot(peak_times, f(peak_times), 'r--', 'linewidth', 1.5);  % Fitted curve in red dashed line
    hold on;

    % Get labels and GOF (Goodness of Fit)
    stim_label = new_out.bipolar_labels{stim};
    resp_label = new_out.bipolar_labels{resp};
    
    % Set the axis limits for consistent annotation placement
    xlim(zoom_times);
    ylim([-1.2 * max(abs(avg)), 1.2 * max(abs(avg))]);  % Dynamic y-limits based on signal amplitude
    
    % Place the annotation at the top-left corner of the plot
    xl = xlim;  % Get the current x-axis limits
    yl = ylim;  % Get the current y-axis limits
    text(xl(1) + 0.02 * (xl(2) - xl(1)), yl(2) - 0.05 * (yl(2) - yl(1)), ...
        sprintf('Stim: %s\nResp: %s\nGOF: %.2f\na: %.2f\nb: %.2f', stim_label, resp_label, gof.adjrsquare, a, b), ...
        'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', 9);

    % Mark the stimulation point with a dashed line at time 0
    plot([0 0], ylim, 'k--');
end