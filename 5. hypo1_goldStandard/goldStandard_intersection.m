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
    ccep_overlap_stims = df{i, 7}{1};
    ccep_overlap_resps = df{i, 8}{1};
    ccep_only_stims = df{i, 9}{1};
    ccep_only_resps = df{i, 10}{1};

    % retrieve the corresponding patient csv out file
    pt_id = df{i, 1}{1};
    sz_id = df{i, 2}{1};
    patient_file = fullfile(thirdOut_dir, [pt_id, '.mat']);
    temp = load(patient_file);
    out = temp.new_out;

    % func to plot 
    num_subplots = size(df{i,6}{1}, 2);
    ccep_assessment(out, num_subplots, pt_id, sz_id, ccep_overlap_stims, ccep_overlap_resps, ccep_only_stims, ccep_only_resps)
    
end