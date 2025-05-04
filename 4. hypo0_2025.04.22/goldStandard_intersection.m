%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
thirdOut_dir = locations.thirdOut_dir;
csv_dir = locations.sz_ccep_overlap;

% load table for ccep spread overlap
df = csv2tbl(csv_dir);

% summary stats of df
num_sz = height(df);
num_pt = length(unique(df{:,1}));

% for each row/seizure in sz_ccep_overlap
for i = 1:num_sz

    % obtain the set of ccep spread channels that had & hadn't intersect
    % with seizure spread channel (i.e. the ccep spread that took the same
    % route as seizure spread and those that had spread to a different
    % region)
    stim_resp.ccep_sig_overlap_stims = df{i, 9}{1};
    stim_resp.ccep_sig_overlap_resps = df{i, 10}{1};
    stim_resp.ccep_sig_only_stims = df{i, 11}{1};
    stim_resp.ccep_sig_only_resps = df{i, 12}{1};
    stim_resp.ccep_nonsig_overlap_stims = df{i, 13}{1};
    stim_resp.ccep_nonsig_overlap_resps = df{i, 14}{1};

    % retrieve the corresponding patient csv out file
    pt_id = df{i, 1}{1};
    sz_id = df{i, 2}{1};
    patient_file = fullfile(thirdOut_dir, [pt_id, '.mat']);
    temp = load(patient_file);
    out = temp.out;

    % func to plot 
    num_subplots = size(df{i,6}{1}, 2) + size(df{i,14}{1}, 2); % all ccep sig + ccep nonsig overlap 
    ccep_assessment(out, num_subplots, pt_id, sz_id, stim_resp)
    
end