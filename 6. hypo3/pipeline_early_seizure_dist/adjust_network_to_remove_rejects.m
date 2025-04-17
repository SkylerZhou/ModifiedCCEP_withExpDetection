%% Aim: for the rejected signal, replace amplitude and latencies of N1 and N2 with NaN. 
function out = adjust_network_to_remove_rejects(out)


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





% commented out/ replacement of mat2csv_hypo2 on 04/16/2025
%{
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

    %% Pick intracranial chs with bipolar signal
    keep_chs = get_chs_to_ignore(out.bipolar_labels); 

    %% Get rejection details arrays
    thresh = out.rejection_details(which_n).thresh;
    which = out.rejection_details(which_n).which;

    sig_avg = out.rejection_details(which_n).reject.sig_avg;
    pre_thresh = out.rejection_details(which_n).reject.pre_thresh;
    at_thresh = out.rejection_details(which_n).reject.at_thresh;
    keep = out.rejection_details(which_n).reject.keep;
    exp = out.rejection_details(which_n).reject.exp;
    ignore_ch = out.rejection_details(which_n).reject.ignore_ch;

    any_reject = sig_avg == 1| pre_thresh == 1 | at_thresh == 1 | exp ==1 | ignore_ch ==1 ;


    %% Restrict to those where rejection_details.keep == 1
    meet_criteria = find(keep==1);
    [row,col] = ind2sub(size(keep),meet_criteria);  % obtain the linear indices of all the keeps 
    % sz comment out on 04/14/2025 with the addition of
    % rejection_details.ignore_ch
    %meet_criteria(keep_chs(row) == false) = []; % filter out the keeps that are recorded with the bad channels(not sure)
    %col(keep_chs(row) == false) = NaN;
    %meet_criteria(keep_chs(col) == false) = [];

    all_stim_idx = unique(row);
    for i = 1:length(all_stim_idx)
        stim_idx = all_stim_idx(i);
        pair_indices = find(row == stim_idx);

        for j = 1:length(pair_indices)
            resp_idx = col(pair_indices(j));
            if ~isnan(resp_idx)
                out.elecs(stim_idx).(adj_n)(resp_idx,:) = out.elecs(stim_idx).(ori_n)(resp_idx,1:2);
            end
        end
    end
end
%}