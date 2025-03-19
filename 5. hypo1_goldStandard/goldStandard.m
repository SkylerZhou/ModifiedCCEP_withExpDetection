% this script aims to plot & humanly-verify the ccep sig and non sig
% spreads

%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
thirdOut_dir = locations.thirdOut_dir;

% load table for ccep sig spread 
spread = readtable('/Users/zhouzican/Documents/CNT/DSOSD/SZ_hypo1_data/goldStandard/sz_ccep_onsetMatches.csv');
num_rows = height(spread);


%% process the network in out so that the out.elecs only display sig vs. nonsig ccep for the new_build_network func  
    
for row=4:4

    % retrieve the patient out file
    patient_id = spread{row, 2}{1};
    patient_file = fullfile(thirdOut_dir, [patient_id, '.mat']);
    temp = load(patient_file);
    out = temp.new_out;


    % retrieve the index of the set of ccep stim channels that have shown seizure actvitiy for this patient
    s = spread{row, 5}{1};
    % convert string into cell array
    ccep_stim = erase(s, {'[', ']'});  
    ccep_stim = strrep(strtrim(strsplit(ccep_stim, ',')), '''', '');  
    ccep_stim = ccep_stim(:);
    % find the indices of these ccep stim chs
    ccep_bipolar_labels = out.bipolar_labels;
    for idx = 1:length(ccep_bipolar_labels)
        if isnumeric(ccep_bipolar_labels{idx}) && isempty(ccep_bipolar_labels{idx})
            ccep_bipolar_labels{idx} = 'NaN';  
        end
    end
    stim_idx = find(ismember(ccep_bipolar_labels, ccep_stim));
    
    
    % use these indices to keep only the corresponding rows in out.elecs 
    out.elecs = out.elecs(stim_idx);
    
    % use the updated out in new_build_network
    save_name = spread{row,3}{1};
    new_build_network(out, 1, save_name)
end

