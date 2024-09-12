overwrite = 1;

%% Updated pipeline to run through all patients in an csv file
locations = cceps_files;
data_folder = locations.data_folder;
results_folder = locations.results_folder;
out_folder = [results_folder,'modified_pipeline/'];

pwfile = locations.pwfile;
login_name = locations.loginname;
script_folder = locations.script_folder;
results_folder = locations.results_folder;

% add paths
addpath(genpath(script_folder));
if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

% Check if output folder exists, if not, create it
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end

%% loop over patients
num_patient = 44;
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));

for n = 44:num_patient
    patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(n));
    out = load(patient_file);
    original_out = out.pt_out;

    %% RUN SINGLE OUT FILE 
    new_out = AA_alternative_filtering(original_out);
    new_out = AA_Running_RejectOrKeep_RW(new_out);
    new_out = AA_new_build_network(new_out); 
    new_out = AA_require_both_Ns(new_out);

    % Save the patient output file
    out_file_name = patient_files(n);
    save(fullfile(out_folder, out_file_name), 'new_out');

    % if there is any error in this step (keeps not sufficient to build the
    % figure), display error and continue
    try
        AA_random_rejections_keeps(new_out);
    catch ME
        fprintf('Error in AA_random_rejections_keeps for patient %s: %s\n', patient_files(n), ME.message);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison between avg and detrend_filt_avg after
% AA_alternative_filtering
num_patient = 19;
row = 133;
col = 50;

for n = 19:num_patient
    patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(n));
    out = load(patient_file);
    out = out.pt_out;

    %% RUN SINGLE OUT FILE 
    out = AA_alternative_filtering(out);
    stim_idx = out.elecs(row).stim_idx;
    if out.elecs(row).N1(col,2) > 0
        ori_n1_idx = out.elecs(row).N1(col,2)+stim_idx;
        ori_n1_x = convert_indices_to_times(ori_n1_idx, out.other.stim.fs, -0.5);
        ori_n1_y = out.elecs(row).detrend_filt_avgs(ori_n1_idx,col);
    else
        ori_n1_idx = 0;
    end
    out = AA_Running_RejectOrKeep_RW(out);
    n1_idx = out.elecs(row).N1(col,2)+stim_idx;
    n1_x = convert_indices_to_times(n1_idx, out.other.stim.fs, -0.5);
    n1_y = out.elecs(row).detrend_filt_avgs(n1_idx,col);
end

% parameters
result = out;
time = convert_indices_to_times(1:1332, result.other.stim.fs, -0.5);
patient = out.name;
stim = out.bipolar_labels{row};
resp = out.bipolar_labels{col};


% plot detrend_filt_avg
tiledlayout(3,1)
nexttile
davg = result.elecs(row).detrend_filt_avgs(:,col);
plot(time, davg, 'black')
xlim([-0.3 0.3])
xline(0,'-r')
detrend_y_lim = ylim;
title_str = sprintf('detrend avg from patient %s with stim at %s and resp at %s', patient, stim, resp);
title(title_str)

% plot n1
hold on
plot(n1_x,n1_y,'x', ...
    'LineWidth',2)
annotation_str = sprintf('N1 (%d)', n1_idx);
text(n1_x, n1_y, annotation_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'red')



% plot avg in the ylim of detrend_filt_avg
nexttile
avg = result.elecs(row).avg(:,col);
plot(time, avg, 'black')
xlim([-0.3 0.3])
xline(0,'-r')
ylim(detrend_y_lim)
title('avg with detrend ylim')

% plot n1
if ori_n1_idx ~= 0
    hold on
    plot(ori_n1_x,ori_n1_y,'x', ...
        'LineWidth',2)
    annotation_str = sprintf('N1 (%d)', ori_n1_idx);
    text(ori_n1_x, ori_n1_y, annotation_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'red')
else
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(2), y_lim(1), 'No N1 detected in original code.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'red')  % Right lower corner
end



% plot avg in regular ylim
nexttile
avg = result.elecs(row).avg(:,col);
plot(time, avg, 'black')
xlim([-0.3 0.3])
xline(0,'-r')
title('avg with regular ylim')

% plot n1
if ori_n1_idx ~= 0
    hold on
    plot(ori_n1_x,ori_n1_y,'x', ...
        'LineWidth',2)
    annotation_str = sprintf('N1 (%d)', ori_n1_idx);
    text(ori_n1_x, ori_n1_y, annotation_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'red')
else
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(2), y_lim(1), 'No N1 detected in original code.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'red')  % Right lower corner
end

