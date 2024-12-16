overwrite = 1;


%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
results_folder = locations.results_folder;
mat_folder = [results_folder,'new_pipeline_keptonly/'];

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));


%% loop over patients
for n = 19:19


    %% load patient out file 
    patient_file = fullfile(mat_folder, patient_files(n));
    temp = load(patient_file);
    out = temp.out;
    %
end