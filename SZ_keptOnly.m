overwrite = 1;

%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
patientNewOut_dir = locations.patientNewOut_dir;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
results_folder = locations.results_folder;
out_folder = [results_folder,'new_pipeline_keptonly/'];

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

% Check if output folder exists, if not, create it
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end
%

%% loop over patients
for n = 1:num_patient

    % load patient out file 
    patient_file = fullfile(patientNewOut_dir, patient_files(n));
    temp = load(patient_file);
    out = temp.new_out;
    
    % adjust amp and lat to remove those that are rejected for n1&n2 
    out = SZ_adjust_network_to_remove_rejects(out); 

    % save the patient output file
    out_file_name = patient_files(n);
    save(fullfile(out_folder, out_file_name), 'out');
end

