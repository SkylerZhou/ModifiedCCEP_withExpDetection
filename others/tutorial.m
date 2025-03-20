%% initial setup 

% step 1. (note that you can skip step 1 if you do not wish to directly run the
% pipeline but only use the patients' files.)
% 1.1. before running the code, get access to the ieeg HUP dataset. download the ieeg.org toolbox to be able to interface with ieeg.org
% 1.2. set up your own directories in cceps_files_example.m. 

% step 2. 
% set up the rest of the directories in cceps_files_examples.m. 
%


%% to open the mat out file for one/multiple patients
locations = cceps_files; % change the name of cceps_files_example.m to cceps_files.m once you have done editing it
data_folder = locations.data_folder;
results_folder = locations.results_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));
%out_folder = [results_folder,'third_pipeline/']; % change third_pipeline/ to your own output directory 
out_folder = [results_folder,'new_pipeline_keptonly/'];

for n = 1:1     % if only want to load the first patient's out file 
    patient_file = fullfile(out_folder, patient_files(n));
    temp = load(patient_file);
    out = temp.out;
end

%{
num_patient = height(ptT);
for n = 1:num_patient     % if only want to loop through all the patients' out file 
    patient_file = fullfile(out_folder, patient_files(n));
    temp = load(patient_file);
    out = temp.out;
end
%}
%


%% to examine the out file (if you are looking at the out file of patient 211)
out.elecs(65).N1(:,1)  % N1 amplitude z-scores of all the responding electrodes to stimulating electrode 65th.
out.elecs(65).n1_adj(:,1)  % Adjusted N1 amplitude z-scores of all the responding electrodes to stimulating electrode 65th. The signals rejected by the CCEP detection algorithm give NaN values. 
out.elecs(77).N1(:,2)  % N1 latency index of all the responding electrodes to stimulating electrode 77th.
out.elecs(77).n1_adj(:,2)  % N1 latency index of all the responding electrodes to stimulating electrode 77th. The signals rejected by the CCEP detection algorithm give NaN values.
% similarly, the above data for N2 can be found in 
% out.elecs(stim_elec_num).N2 
% out.elecs(stim_elec_num).n2_adj
%


%% to visualize the adjancency matrix
show_network(out,1,1,0) 
new_build_network(out,1)
%


%% to view the reasons for ccep rejections 
out.rejection_details(1).reject
% sig_avg, pre_thresh, at_thresh were reasons of rejections defined by Erin
% no_n1, no_both, empty were set by Rudy
% exp was set by skyler


%% to run Erin's pipeline (the first version) 
% download Erin's pipeline from https://github.com/erinconrad/CCEPS/tree/main
% set up the ccep_files.m in the folder holding Erin's script
% run the pipeline using do_run/updated_pipeline/do_all_pts_v2.m
% the result from this pipeline should be stored at your locations.results_folder


%% to run Skyler's pipeline (the third version, which had already incorporated Rudy's (second) pipeline in 1.compare_pipelines so there is no need to use the output from Erin's pipeline and loop them through Rudy's pipeline  
% use the outputs from Erin's pipeline in result_folder to run the third
% pipeline using this portion of code in SZ_runNew_compareOriNew.m
overwrite = 1;
set(0, 'DefaultFigureVisible', 'off'); % prevent figures from poping up

% prep
% dir to input files and scripts
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));
script_folder = locations.script_folder;

% to acess the 1st version of patients' output
firstOut_dir = locations.firstOut_dir;

% to store the output of the 3rd version
thirdOut_dir = locations.thirdOut_dir;

% add ieeg paths
pwfile = locations.pwfile;
login_name = locations.loginname;

addpath(genpath(script_folder));
if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

% loop over patients
start_patient = 1;
num_patient = 1;
for n = start_patient:num_patient
    patient_file = fullfile(firstOut_dir, patient_files(n));
    out = load(patient_file);
    ori_out = out.pt_out;

    % run single out file 
    new_out = RW_alternative_filtering(ori_out); 
    new_out = RW_Running_RejectOrKeep(new_out); 
    new_out = RW_new_build_network(new_out); 
    new_out = RW_require_both_Ns(new_out);
    % SZ_adjust_network_to_remove_rejects replace zeros and/or n1&n2 amplitudes with NaN if they are deemed as non-sig ccep signal by rudy's & skyler's versions of code
    % will only be executed on two newly constructed fields (out.elecs.n1_adj and
    % out.elecs.n2_adj), meaning that running this code will not affect
    % out.elecs.N1 and out.elecs.N2
    new_out = SZ_adjust_network_to_remove_rejects(new_out);  

    % save the patient output file
    out_file_name = patient_files(n);
    save(fullfile(thirdOut_dir, out_file_name), 'new_out');

    % RW_random_rejections_keeps display a random set (n=25) of significant and
    % non-significant ccep signals so that you can verify if the kept ones
    % are true positives and if the rejected ones are true negatives. 
    % I set set(0, 'DefaultFigureVisible', 'off'); to prevent figures from
    % poping up. You can change this to on to view the output from RW_random_rejections_keeps
    try
        RW_random_rejections_keeps(new_out);
    catch ME
        fprintf('Error in random_rejections_keeps for patient %s: %s\n', patient_files(n), ME.message);
    end
end
%