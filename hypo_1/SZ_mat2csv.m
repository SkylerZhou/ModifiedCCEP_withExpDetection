overwrite = 1;


%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
results_folder = locations.results_folder;
mat_folder = [results_folder,'new_pipeline_keptonly/'];
csv_folder = [results_folder,'new_pipeline_keptonly_csv/'];

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

% Check if the csv output folder exists, if not, create it
if ~exist(csv_folder, 'dir')
    mkdir(csv_folder);
end
%


%% loop over patients
for n = 52:52


    %% load patient out file 
    patient_file = fullfile(mat_folder, patient_files(n));
    temp = load(patient_file);
    out = temp.out;
    %

    %% extract N1, stim, and response channels info 
    num_elecs = size(out.elecs, 2);
    n1_all = [];
    stim_all = [];
    resp_all = [];
    
    for ich = 1:num_elecs
        % if the current elecs is not a stimulating elecs, skip it
        if isempty(out.elecs(ich).arts), continue; end

        % if exist distance info, extract n1 amp and lat, and distance
        if isfield(out.other, 'elecs_dist') 
            % cp the distance info as the thrid column to n1_adj
            out.elecs(ich).n1_adj(:,3) = out.other.elecs_dist(:,ich);
            % extract n1 amp and lat, and distance
            n1 = out.elecs(ich).n1_adj;
            n1_all = [n1_all; n1];
        else
            % else extract only n1 amp and lat
            n1 = out.elecs(ich).n1_adj;
            n1_all = [n1_all; n1];

        end

        % extract stim and response labels 
        stim = repmat(out.bipolar_labels(ich), size(out.elecs(ich).n1_adj,1), 1);
        stim_all = [stim_all; stim];
        resp = out.bipolar_labels(1: size(out.elecs(ich).n1_adj,1));
        resp_all = [resp_all; resp];
    end



    %% create cell arrays containing patient, stim_channel, response_channel, n1_amp, n1_lat, elecs_dist
    num_rows = size(n1_all, 1);
    num_cols = 6; % patient, stim_channel, response_channel, elecs_dist, n1_amp, n1_lat
    new_out = cell(num_rows, num_cols);
    
    % add patient id to col 1
    patient_col = repmat({out.name}, num_rows, 1);
    new_out(:,1) = patient_col;
    
    % add stim_channel and resp_channel to col 2 and 3
    new_out(:,2) = stim_all;
    new_out(:,3) = resp_all;

    % append elecs_dist, n1_amp, n1_lat to col 4 to 6
    n1_all = num2cell(n1_all);
    if size(n1_all, 2) == 3
        % If distance information is available
        new_out(:,5) = n1_all(:,1); % n1 amp
        new_out(:,6) = n1_all(:,2); % n1 lat
        new_out(:,4) = n1_all(:,3); % elecs_dist
    else
        % If distance information is not available
        new_out(:,5) = n1_all(:,1); % n1 amp
        new_out(:,6) = n1_all(:,2); % n1 lat
        new_out(:,4) = num2cell(NaN(num_rows,1)); % Placeholder for elecs_dist
    end
    %


    %% save the new output cell array file as csv
    varNames = {'PatientID', 'StimChannel', 'RespChannel', 'ElecsDist', 'N1Amp', 'N1Lat'};
    T = cell2table(new_out, 'VariableNames', varNames);

    % write table to CSV
    output_file = fullfile(csv_folder, [out.name, '.csv']);
    writetable(T, output_file);
    fprintf('Saved csv for %s', out.name);
    %
end

