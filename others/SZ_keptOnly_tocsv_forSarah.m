overwrite = 1;

%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
results_folder = locations.results_folder;
file_folder = [results_folder,'new_pipeline_keptonly/'];
csv_out_folder = [results_folder,'new_pipeline_keptonly_csv/'];

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

% Check if the csv output folder exists, if not, create it
if ~exist(csv_out_folder, 'dir')
    mkdir(csv_out_folder);
end
%

%% loop over patients
%for n = 1:num_patient
for n = 2:2   

    % load patient out file 
    patient_file = fullfile(file_folder, patient_files(n));
    temp = load(patient_file);
    out = temp.out;

    for ich = 1:num_elecs
        % if the current elecs is not a stimulating elecs, skip it as it
        % contains no data to be altered.
        if isempty(out.elecs(ich).arts), continue; end
            
        % extract n1 
        n1 = out.elecs(ich).N1(:,1:2);
    end
    
    % create empty matrix with size = (number of stim * number of resp) rows * 6 columns
    num_cols = 6; % patient, stim_channel, response_channel, distance_between_channels, N1_amp, N1_lat
    num_stims = size(out.other.stim_elecs,1);
    num_resp = size(out.bipolar_ch_pair,1);

    
    % if has elecs distance info
    % out.other.elecs_dist

    % save the patient output file
    % out_file_name = patient_files(n);
    % save(fullfile(out_folder, out_file_name), 'out');
end
