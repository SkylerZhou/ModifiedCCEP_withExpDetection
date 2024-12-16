%% initial setup 

% step 1. (note that you can skip step 1 if you do not wish to directly run the
% pipeline but only use the patients' files.)
% 1.1. before running the code, get access to the ieeg HUP dataset. download the ieeg.org toolbox to be able to interface with ieeg.org
% 1.2. set up your own directories in cceps_files.m. Specifically, 
% 1.2.1. locations.pwfile: stores the directory to your _ieeglogin.bin pwd
% 1.2.2. locations.loginname: add your login name to ieeg.org

% step 2. 
% 2.1. set up your own directories in cceps_files.m. Specifically,
% 2.1.1. locations.data_folder: stores a pt_mat file that you can download in https://docs.google.com/spreadsheets/d/1qXCDxgi9CC3ZN--OxNMukBWiGDBvhyxok_9VuzP0CD4/edit?usp=sharing
% 2.1.2. locations.result_folder: stores the directory where you want to store the output in
% 2.1.3. locations.script_folder: stores your own scripts such as tutorial.m
% or if you want to use any other functions that Erin, Rudy, or I wrote
% such as show_network()
% 2.1.4. you can ignore/delete the rest such as locations.pwfile and locations.patientNewOut_dir.
%


%% to open the mat out file for one/multiple patients
locations = cceps_files;
data_folder = locations.data_folder;
results_folder = locations.results_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));
out_folder = [results_folder,'new_pipeline_keptonly/']; % change the new_pipeline_keptonly to your own directory name

for n = 1:1     % if only want to load the first patient's out file 
    patient_file = fullfile(out_folder, patient_files(n));
    temp = load(patient_file);
    out = temp.out;
end

%{
num_patient = height(ptT);
for n = 1:num_patient     % if only want to load all the patients' out file 
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
% similarly, the above data for N2 can be found in out.elecs(stim_elec_num).N2 and
% out.elecs(stim_elec_num).n2_adj
%


%% to visualize the adjancency matrix
show_network(out,1,1,0)
%