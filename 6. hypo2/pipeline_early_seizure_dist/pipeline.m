% this script aims to:
% 1. add_elecs_distance - add the elec distance info to each pair of stim-resp pair
% 2. add_roi - add the region info 
% 3. adjust_network_to_remove_rejects - adjust the n1_adj and n2_adj to formulate a new csv
% output. The peak amp of n1_adj and n2_adj in this version of should only == NaN if its corresponding
% sig_avg == 1 | pre_thresh == 1 | exp == 1 | ignore_ch ==1 . 
% 4. mat2csv - generate csv for python processing

overwrite = 1;

%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
thirdOut_dir = locations.thirdOut_dir;
csv_folder = [thirdOut_dir,'third_pipeline_hypo2/early_seizure_dist/'];

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

    %% try add_elecs_distance function; check if there are corresponding
    % coordinate files to build distance matrix for this patient 
    try 
        out = calculate_bipolar_midpoint(out);
        out = add_bipolar_dist(out);
    catch ME 
        fprintf('%s\n', ME.message);
    end
    
    %% adjust amp and lat to remove those that are rejected for n1&n2 
    %out = adjust_network_to_remove_rejects(out); 

    % save the patient output file
    %out_file_name = patient_files(n);
    %save(fullfile(thirdOut_dir, out_file_name), 'out');



    %% if patient's ccep output performance is bad (through visual examination), exclude these patients
    pt_id = out.name;
    bad_data = ["HUP213", "HUP214", "HUP216", "HUP256", "HUP264", "HUP266", "HUP272", "HUP273"];
    if any(strcmp(pt_id, bad_data))
        fprintf('Skipping %s due to poor CCEP performance.\n', pt_id);
        continue;  % skip to next patient
    end


    %% convert to csv for processsing in python
    [T, csv_name] = mat2csv(out, csv_folder);
    

    % save the csv 
    writetable(T, csv_name);
    fprintf('Saved csv for %s', out.name);
end