%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comparison between the detrend_filt_avg and the detrend_filt_avg after interpolation
patient_idx = 19;
row = 133;
col = 50;
%


%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
patientOriOut_dir = locations.patientOriOut_dir;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));

% retrieve the ori_out
patient_file = fullfile(patientOriOut_dir, patient_files(patient_idx));
temp = load(patient_file);
out = temp.pt_out;
% 



%% run ori_out through alternative_filtering to obtain the detrend_filt_avg
% that has not undergone interpolation yet
out = RW_alternative_filtering(out);
stim_idx = out.elecs(row).stim_idx;

if out.elecs(row).N1(col,2) > 0
    ori_n1_idx = out.elecs(row).N1(col,2)+stim_idx;
    ori_n1_x = convert_indices_to_times(ori_n1_idx, out.other.stim.fs, -0.5);
    ori_n1_y = out.elecs(row).detrend_filt_avgs(ori_n1_idx,col);
else
    ori_n1_idx = 0;
end


% run the output through running_rejectorkeep to obtain the
% detrend_filt_avg after interpolation
out = RW_Running_RejectOrKeep(out);
new_n1_idx = out.elecs(row).N1(col,2)+stim_idx;
new_n1_x = convert_indices_to_times(new_n1_idx, out.other.stim.fs, -0.5);
new_n1_y = out.elecs(row).detrend_filt_avgs(new_n1_idx,col);
%



%% plot
% parameters
result = out;
patient = out.name;
time = convert_indices_to_times(1:1332, result.other.stim.fs, -0.5);
stim = out.bipolar_labels{row};
resp = out.bipolar_labels{col};


%% plot detrend_filt_avg after interpolation
tiledlayout(3,1)
nexttile
davg_int = result.elecs(row).detrend_filt_avgs(:,col);
plot(time, davg_int, 'black')
xlim([-0.3 0.3])
xline(0,'-r')
davg_int_ylim = ylim;
title_str = sprintf('Detrend avg after interpolation from patient %s with stim at %s and resp at %s', patient, stim, resp);
title(title_str)

% plot n1
hold on
plot(new_n1_x,new_n1_y,'x', ...
    'LineWidth',2)
annotation_str = sprintf('N1 (%d)', new_n1_idx);
text(new_n1_x, new_n1_y, annotation_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'red')



%% plot detrend_filt_avg using the scale of davg_int_ylim
nexttile
davg = result.elecs(row).avg(:,col);
plot(time, davg, 'black')
xlim([-0.3 0.3])
xline(0,'-r')
ylim(davg_int_ylim)
title_str = sprintf('Detrend avg without interpolation from patient %s with stim at %s and resp at %s', patient, stim, resp);
title(title_str)


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



%% plot detrend_filt_avg in regular ylim
nexttile
davg = result.elecs(row).avg(:,col);
plot(time, davg, 'black')
xlim([-0.3 0.3])
xline(0,'-r')
title_str = sprintf('Detrend avg without interpolation from patient %s with stim at %s and resp at %s', patient, stim, resp);
title(title_str)

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

