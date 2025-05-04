% this script aims to:
% 1. add the elec distance info to each pair of stim-resp pair 
% 2. adjust the n1_adj and n2_adj to formulate a new csv
% output. The peak amp of n1_adj and n2_adj in this version of should only == NaN if its corresponding
% sig_avg == 1 | pre_thresh == 1 | exp == 1 | ignore_ch ==1 . 
% 3. check if the ccep performance for this patient in the detector is bad,
% if bad, do not save it as csv to aviod downstream hypothesis testing.
% 3. generate csv for python processing

overwrite = 1;

%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
thirdOut_dir = locations.thirdOut_dir;
hypo2_dir = [thirdOut_dir,'third_pipeline_hypo2/elecDist/'];
csv_folder = [hypo2_dir,'csv/'];

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

    % Remove old field if it exists
    if isfield(out, 'other') && isfield(out.other, 'elecs_dist')
        out.other = rmfield(out.other, 'elecs_dist');
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
    out = adjust_network_to_remove_rejects(out); 

    % save the patient output file
    out_file_name = patient_files(n);
    save(fullfile(hypo2_dir, out_file_name), 'out');


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




% commented out on 04/16/2025 to to make add_elecs_distance,
% adjust_network_to_remove_reject, and mat2csv functions - break the long
% func to make it clearer. 
%{
%% old code 

% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);
thirdOut_dir = locations.thirdOut_dir;
csv_folder = [thirdOut_dir,'third_pipeline_csv_hypo2/'];

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

% Check if the csv output folder exists, if not, create it
if ~exist(csv_folder, 'dir')
    mkdir(csv_folder);
end
%

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

    %% Pre-allocate n1_adj and n2_adj matrix with nan
    num_elecs = size(out.elecs, 2); % get num of electrodes for this patient 
    for row = 1:num_elecs
        if size(out.elecs(row).avg, 1) >= 1
            out.elecs(row).n1_adj = nan(num_elecs, 2); 
            out.elecs(row).n2_adj = nan(num_elecs, 2); 
        end
    end

    ori = {'N1', 'N2'};
    adj = {'n1_adj','n2_adj'};

    for j = 1:2
        if j == 1
            which_n = 1; % interested in n1
            ori_n = ori{j};
            adj_n = adj{j};
        else
            which_n = 2; % n2
            ori_n = ori{j};
            adj_n = adj{j};
        end

        %% Get rejection details arrays
        %thresh = out.rejection_details(which_n).thresh;
        %which = out.rejection_details(which_n).which;
        %at_thresh = out.rejection_details(which_n).reject.at_thresh;
        %keep = out.rejection_details(which_n).reject.keep;
        sig_avg = out.rejection_details(which_n).reject.sig_avg;
        pre_thresh = out.rejection_details(which_n).reject.pre_thresh;
        exp = out.rejection_details(which_n).reject.exp;
        ignore_ch = out.rejection_details(which_n).reject.ignore_ch;
        good = (sig_avg ~= 1) & (pre_thresh ~= 1) & (exp ~= 1) & (ignore_ch ~= 1);

        %% leave stim-resp peak amp to nan if not in good 
        good_indices = find(good==1);
        [row,col] = ind2sub(size(good),good_indices);  % obtain the linear indices of all the keeps     
        all_stim_idx = unique(row);

        for i = 1:length(all_stim_idx)
            stim_idx = all_stim_idx(i);
           
            if size(out.elecs(stim_idx).arts,1) ~= 0 % check if the ch is actually stimulating
                pair_indices = find(row == stim_idx); % find the resp for those stim ch
    
                for j = 1:length(pair_indices)
                    resp_idx = col(pair_indices(j));
                    out.elecs(stim_idx).(adj_n)(resp_idx,:) = out.elecs(stim_idx).(ori_n)(resp_idx,1:2);          
                end
            end
        end
   
    end



    %% covert to csv 
    % extract N1, stim, and response channels info 
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
%}