% Replace non-string entries with empty strings
search_out = ori_out.bipolar_labels;
is_nonstring = cellfun(@(x) ~ischar(x), search_out);
search_out(is_nonstring) = {''};

% set stim and resp
search_stim = 'LB6-LB7';
search_resp = 'LB2-LB3';

% Use strcmp to compare each element in the cell array with the search term
stim_matches = strcmp(search_out, search_stim);
resp_matches = strcmp(search_out, search_resp);

% Find the indices of the matching entries
[stim_row, stim_col] = find(stim_matches);
[resp_row, resp_col] = find(resp_matches);

% Display the result
if isempty(stim_row)
    disp('Search term not found.');
else
    fprintf('Stim: %d\nResp: %d', stim_row, resp_row);
end



%%%%%%
times = new_out.elecs(stim_row).times;
davg = new_out.elecs(stim_row).detrend_filt_avgs(:,resp_row);
eeg_times = convert_indices_to_times(1:length(davg),new_out.other.stim.fs,times(1));

plot(eeg_times,davg,'k','linewidth',2);
zoom_times = [-300e-3 300e-3];
xlim([zoom_times(1) zoom_times(2)]);