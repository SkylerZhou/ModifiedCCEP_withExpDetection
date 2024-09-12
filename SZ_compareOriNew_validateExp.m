% initialize arrays to store the total number of keeps for erin's and
% rudy's algo
patients = cell(52, 1);
ori_keep = zeros(52, 1);
new_keep = zeros(52, 1);
total_response = zeros(52, 1);

%% loop over patients
num_patient = 23;
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));

for n = 23:num_patient
    patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(n));
    out = load(patient_file);
    original_out = out.pt_out;

    %% RUN SINGLE OUT FILE 
    new_out = RW_alternative_filtering(original_out);
    new_out = RW_Running_RejectOrKeep(new_out);
    new_out = RW_new_build_network(new_out); 
    new_out = RW_require_both_Ns(new_out);

    % if there is any error in this step (keeps not sufficient to build the
    % figure), display error and continue
    try
        RW_random_rejections_keeps(new_out);
    catch ME
        fprintf('Error in random_rejections_keeps for patient %s: %s\n', patient_files(n), ME.message);
    end

    
    %% yield the total number of keeps and rejects for the erin's original_out
    which_n = 1;
    keep_ori = original_out.rejection_details(which_n).reject.keep;
    sig_avg_ori = original_out.rejection_details(which_n).reject.sig_avg;
    pre_thresh_ori = original_out.rejection_details(which_n).reject.pre_thresh;
    at_thresh_ori = original_out.rejection_details(which_n).reject.at_thresh;
    
    any_reject_ori = sig_avg_ori == 1 | pre_thresh_ori == 1 | at_thresh_ori == 1;

    original_nkeep = sum(keep_ori(:) == 1);
    original_nreject = sum(any_reject_ori(:) == 1);

    % store the ori keeps in array 
    patients{n} = char(ptT.HUPID{n});
    ori_keep(n) = original_nkeep;


    %% yield the total number of keeps and rejects for the rudy's new_out
    keep_new = new_out.rejection_details(which_n).reject.keep;
    sig_avg_new = new_out.rejection_details(which_n).reject.sig_avg;
    pre_thresh_new = new_out.rejection_details(which_n).reject.pre_thresh;
    at_thresh_new = new_out.rejection_details(which_n).reject.at_thresh;
    no_both = new_out.rejection_details(which_n).reject.no_both;
    exp = new_out.rejection_details(which_n).reject.exp;

    any_reject_new = sig_avg_new == 1 | pre_thresh_new == 1 | at_thresh_new == 1 | no_both == 1 | exp == 1;

    new_nkeep = sum(keep_new(:) == 1);
    new_nreject = sum(any_reject_new(:) == 1);
    
    %if original_nkeep + original_nreject ~= new_nkeep + new_nreject
    %    continue
    %end 
    % HUP231 and HUP239 total number do not match
    % store the new keeps in array
    new_keep(n) = new_nkeep;
    total_response(n) = original_nkeep + original_nreject;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% crate a table from the arrays to compare the total keeps of Erin's and Rudy's algorithms 
ori_keep_perct = rdivide(ori_keep, total_response) * 100;
new_keep_perct = rdivide(new_keep, total_response) * 100;
compare_table = table(patients, total_response, ori_keep, new_keep,ori_keep_perct, new_keep_perct,...
    'VariableNames',{'HUPID', 'Total_Responses', 'Original_Total_Keep', 'New_Total_Keep', 'Original_Perct(%)', 'New_Perct(%)'});
% save locally
writetable(compare_table, 'Compare_Total_Keep.xlsx')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randomly select 20 out of all 52 patients' data for validation
random_select = 20;
filter_table = compare_table(compare_table.New_Total_Keep >= 25 & ...
                               ~strcmp(compare_table.HUPID, 'HUP216'), :);
random_indices = randperm(height(filter_table), random_select);
select_patients = filter_table.HUPID(random_indices);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for those rejected by exp, randomly select and plot 25 figures
% to see how good curve fitting gof is in rejecting true artifacts 
% Parameters
number_to_plot = 25;
zoom_times = [-300e-3 300e-3];

% find row and col for waves that were only rejected due to 'exp' in
% 'new_out' but not rejected in 'original_out'
exp_reject = (exp == 1); 
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
    peak_end_time = 0.3;
    stim_idx = new_out.elecs(stim).stim_idx;
    peak_start_index = new_out.elecs(stim).N1(resp,2) + stim_idx;  % the index of 0 sec is 512
    peak_end_index = stim_idx + floor(new_out.other.stim.fs * peak_end_time);                  
    peak_times = convert_indices_to_times(peak_start_index:peak_end_index, new_out.other.stim.fs, times(1));
    after_peak_avg = new_out.elecs(stim).detrend_filt_avgs(peak_start_index:peak_end_index, resp);                 
    
    % Examine goodness of fit using exponential fit
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'Lower', [-2000, -50], ...
                        'Upper', [2000, 0]);
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