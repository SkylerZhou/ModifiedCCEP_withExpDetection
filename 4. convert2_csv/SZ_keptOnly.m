overwrite = 1;

%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
results_folder = locations.results_folder;
thirdOut_dir = locations.thirdOut_dir;
out_folder = [results_folder,'new_pipeline_keptonly/'];

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));


%% loop over patients
for n = 1:num_patient

    % load patient out file 
    patient_file = fullfile(thirdOut_dir, patient_files(n));
    temp = load(patient_file);
        
    if isfield(temp, 'new_out')
        out = temp.new_out;
    elseif isfield(temp, 'out')
        out = temp.out;
    end

    % try add_elecs_distance function; check if there are corresponding
    % coordinate files to build distance matrix for this patient 
    try 
        out = SZ_add_elecs_distance(out);
    catch ME 
        fprintf('%s\n', ME.message);
    end
    
    % adjust amp and lat to remove those that are rejected for n1&n2 
    out = SZ_adjust_network_to_remove_rejects(out); 

    % save the patient output file
    out_file_name = patient_files(n);
    save(fullfile(thirdOut_dir, out_file_name), 'out');
end

