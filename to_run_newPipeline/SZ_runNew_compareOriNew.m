overwrite = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Updated pipeline to run through all patients in an csv file

%% prep
locations = cceps_files;
data_folder = locations.data_folder;
results_folder = locations.results_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));
out_folder = [results_folder,'new_pipeline/'];
patientNewOut_dir = locations.patientNewOut_dir;
patientOriout_dir = locations.patientOriOut_dir;

% add ieeg paths
pwfile = locations.pwfile;
login_name = locations.loginname;
script_folder = locations.script_folder;
results_folder = locations.results_folder;

addpath(genpath(script_folder));
if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

% Check if output folder exists, if not, create it
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end
%


%% loop over patients
start_patient = 55;
num_patient = 56;

for n = start_patient:num_patient
    patient_file = fullfile(patientOriout_dir, patient_files(n));
    out = load(patient_file);
    ori_out = out.pt_out;

    % run single out file 
    new_out = RW_alternative_filtering(ori_out);
    new_out = RW_Running_RejectOrKeep(new_out);
    new_out = RW_new_build_network(new_out); 
    new_out = RW_require_both_Ns(new_out);

    % save the patient output file
    out_file_name = patient_files(n);
    save(fullfile(out_folder, out_file_name), 'new_out');

    % if there is any error in this step (keeps not sufficient to build the
    % figure), display error and continue
    try
        RW_random_rejections_keeps(new_out);
    catch ME
        fprintf('Error in random_rejections_keeps for patient %s: %s\n', patient_files(n), ME.message);
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% total number of keeps remained after running original and new algorithms
% initialize arrays to store 
arr_patients = cell(58, 1);
arr_sum_keep_ori = zeros(58, 1);
arr_sum_keep_new = zeros(58, 1);
arr_total = zeros(58, 1);


%% loop over to retrieve patient's data
start_patient = 1;
num_patient = 58;
which_n = 1;

for n = start_patient:num_patient

    % load the mat files outputed by running the two versions of codes
    ori_patient_file = fullfile(patientOriOut_dir, patient_files(n));
    temp = load(ori_patient_file);
    ori_out = temp.pt_out;

    new_patient_file = fullfile(patientNewOut_dir, patient_files(n));
    temp = load(new_patient_file);
    new_out = temp.new_out;
 

    % obtain the total number of keeps from ori 
    keep_ori = ori_out.rejection_details(which_n).reject.keep;
    sig_avg_ori = ori_out.rejection_details(which_n).reject.sig_avg;
    pre_thresh_ori = ori_out.rejection_details(which_n).reject.pre_thresh;
    at_thresh_ori = ori_out.rejection_details(which_n).reject.at_thresh;
    
    any_reject_ori = sig_avg_ori == 1 | pre_thresh_ori == 1 | at_thresh_ori == 1;

    sum_keep_ori = sum(keep_ori(:) == 1);
    sum_reject_ori = sum(any_reject_ori(:) == 1);


    % obtain the total number of keeps from new 
    keep_new = new_out.rejection_details(which_n).reject.keep;   
    sig_avg_new = new_out.rejection_details(which_n).reject.sig_avg;
    pre_thresh_new = new_out.rejection_details(which_n).reject.pre_thresh;
    at_thresh_new = new_out.rejection_details(which_n).reject.at_thresh;
    no_both_new = new_out.rejection_details(which_n).reject.no_both;
    exp_new = new_out.rejection_details(which_n).reject.exp;
    
    any_reject_new = sig_avg_new == 1 | pre_thresh_new == 1 | at_thresh_new == 1 | no_both_new == 1 | exp_new == 1;

    sum_keep_new = sum(keep_new(:) == 1);
    sum_reject_new = sum(any_reject_new(:) == 1);
    

    %if original_nkeep + original_nreject ~= new_nkeep + new_nreject
    %    continue
    %end 
    % HUP231 and HUP239 total number do not match
    

    % store the keeps in array
    arr_patients{n} = char(ptT.HUPID{n});
    arr_sum_keep_ori(n) = sum_keep_ori;
    arr_sum_keep_new(n) = sum_keep_new;
    arr_total(n) = sum_keep_ori + sum_reject_ori;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% crate a table from the arrays to compare the total keeps of two versions of algorithms 
ori_keep_perct = rdivide(arr_sum_keep_ori, arr_total) * 100;
new_keep_perct = rdivide(arr_sum_keep_new, arr_total) * 100;
compare_table = table(arr_patients, arr_total, arr_sum_keep_ori, arr_sum_keep_new,ori_keep_perct, new_keep_perct,...
    'VariableNames',{'HUPID', 'Total_Responses', 'Original_Total_Keep', 'New_Total_Keep', 'Original_Perct(%)', 'New_Perct(%)'});
% save locally
writetable(compare_table, 'Compare_Total_Keep.xlsx')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% randomly select 20 out of all 52 patients' data for validation
random_select = 20;
% 214, 216, 256, 266, 272, 273 tend to have poor performances 
% 211, 213, 214, 229 does not yield enough keep validations in Rudy's
% version of code 
filter_table = compare_table(compare_table.New_Total_Keep >= 25 & ...
                               ~strcmp(compare_table.HUPID, 'HUP211') &...
                               ~strcmp(compare_table.HUPID, 'HUP213') &...
                               ~strcmp(compare_table.HUPID, 'HUP214') &...
                               ~strcmp(compare_table.HUPID, 'HUP216') &...
                               ~strcmp(compare_table.HUPID, 'HUP229') &...
                               ~strcmp(compare_table.HUPID, 'HUP256') &...
                               ~strcmp(compare_table.HUPID, 'HUP266') &...
                               ~strcmp(compare_table.HUPID, 'HUP272') &...
                               ~strcmp(compare_table.HUPID, 'HUP273'), :);
random_indices = randperm(height(filter_table), random_select);
select_patients = filter_table.HUPID(random_indices);
