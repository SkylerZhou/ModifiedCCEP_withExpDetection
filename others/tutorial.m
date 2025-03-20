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