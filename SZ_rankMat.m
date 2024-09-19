%% This script aims to construct a connectivity matrix of size (num_elecs * num elecs) for each patient
% 1). Rows are the stimulating electrodes; Columns are the responding
% electrodes.
% 2). For each stimulating electrode, the responding CCEPs will be ranked
% according to their latencies. 
% 3). The ranking will be NaN if no latencies are recorded (or if they are
% rejected) 

% specify directories
num_patient = 32;
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));


% loop over patients
for n = 32:num_patient

    % load patient data 
    new_patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(n));
    temp = load(new_patient_file);
    out = temp.new_out;

    % build matrix to store the rankings 
    num_elecs = size(out.elecs, 2);
    n1_rank_mat = NaN(num_elecs, num_elecs);
    n2_rank_mat = NaN(num_elecs, num_elecs);

    % extract latencies of all responding electrode to a single stimulating
    % electrode
    idx_elecs = 1;
    for ich = 1:num_elecs
        % if the current ich-th elecs is not stimulating, skip it and leave
        % its row in rank_mat all NaN.
        if isempty(out.elecs(ich).arts), continue; end

        %% for N1
        n1_lat = out.elecs(ich).N1(:,2);
        % process to filter out the 0s and the ones rejected 
        % rank according to their latencies
        % store as the idx_elecs row in n1_rank_mat 


    end

end



%% compare two rankMat using Brain Connectivity Toolbox 