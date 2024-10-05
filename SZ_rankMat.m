%% This script aims to construct a connectivity matrix of size (num_elecs * num elecs) for each patient
% 1). Rows are the stimulating electrodes; Columns are the responding
% electrodes.
% 2). For each stimulating electrode, the responding CCEPs will be ranked
% according to their latencies. 
% 3). The ranking will be NaN if no latencies are recorded (or if they are
% rejected) 

%% prep
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
%


%% loop over patients
for n = 18:18

    % load patient out file  
    new_patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(n));
    temp = load(new_patient_file);
    out = temp.new_out;
    % adjust amp and lat for n1&n2 to NaN if they are rejected 
    out = SZ_adjust_network_to_remove_rejects(out); 

    % build matrix to store the rankings 
    num_elecs = size(out.elecs, 2);
    n1_rank_mat = NaN(num_elecs, num_elecs);
    n2_rank_mat = NaN(num_elecs, num_elecs);
    idx_elecs = 1;
    % 

    % loop over all elecs of this patient 
    % extract latencies of all responding elecs of a stim electrode
    for ich = 1:num_elecs
        % if the current ich-th elecs is not stimulating, skip it and leave
        % its row in rank_mat all NaN.
        if isempty(out.elecs(ich).arts), continue; end

        %% for N1
        % rank according to their latencies
        n1_lat = out.elecs(ich).n1_adj(:,2);
        rank_n1_lat = tiedrank(n1_lat, 'nanflag', 'omitnan');

        % store as the idx_elecs row in n1_rank_mat
        n1_rank_mat(ich,:) = rank_n1_lat;


    end

end


%% compare two rankMat using Brain Connectivity Toolbox 