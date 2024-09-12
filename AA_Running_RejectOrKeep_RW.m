% Further processes the data by running a rejection or keep analysis. 
% It examines the filtered data for artifact signals or low-quality
% data points and made a decision on keep or reject.
function out = AA_Running_RejectOrKeep_RW(out)

%% Parameters
idx_before_stim = 30;
n1_time = [11e-3 50e-3];
loose_n1_time = [11e-3 50e-3];
n2_time = [50e-3 300e-3];
stim_time = [-5e-3 10e-3];
tight_stim_time = [-5e-3 10e-3];
stim_val_thresh = 1e3;
rel_thresh = 3;
fs = out.other.stim.fs;
max_crossings = 3;

n1_idx = floor(n1_time*fs);
loose_n1_idx = floor(loose_n1_time*fs);
n2_idx = floor(n2_time*fs);
stim_indices = floor(stim_time*fs);
tight_stim_indices = floor(tight_stim_time*fs);


%{
if size(out.elecs(5).stim_idx,1)~=0
    stim_start = out.elecs(5).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)~=0
    stim_start = out.elecs(7).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)~=0
    stim_start = out.elecs(11).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)~=0
    stim_start = out.elecs(85).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)~=0
    stim_start = out.elecs(18).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)~=0
    stim_start = out.elecs(9).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)~=0
    stim_start = out.elecs(10).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)~=0
    stim_start = out.elecs(23).stim_idx-3;
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)==0 && size(out.elecs(3).stim_idx,1)~=0
    stim_start = out.elecs(3).stim_idx-3;
end
if size(out.chLabels,1)>=185 
    if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)==0 && size(out.elecs(3).stim_idx,1)==0 && size(out.elecs(185),1)~=0
        stim_start = out.elecs(185).stim_idx-3;
    end
end
if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)==0 && size(out.elecs(3).stim_idx,1)==0 && size(out.elecs(185).stim_idx,1)==0 && size(out.elecs(13),1)~=0
    stim_start = out.elecs(13).stim_idx;
end
%}


% Loop over elecs
n = size(out.chLabels,1);
%% 

for ich = 1:n
    if isempty(out.elecs(ich).arts), continue; end
    
    n1 = zeros(size(out.elecs(ich).detrend_filt_avgs,2),4);
    n2 = zeros(size(out.elecs(ich).detrend_filt_avgs,2),4);
    
    stim_idx = out.elecs(ich).stim_idx;
    stim_start = stim_idx;
    
    % redefine n1 and n2 relative to beginning of eeg
    temp_n1_idx = n1_idx + stim_idx - 1;
    temp_loose_n1_idx = loose_n1_idx + stim_idx - 1;
    temp_n2_idx = n2_idx + stim_idx - 1;
    temp_stim_idx = stim_indices + stim_idx - 1;
    temp_tight_stim = tight_stim_indices + stim_idx-1;
    
    
    % Loop over channels within this elec
    for jch = 1:size(out.elecs(ich).detrend_filt_avgs,2)

        % Get the eeg
        eeg = out.elecs(ich).detrend_filt_avgs(:,jch);
        stim_channel = ich;
        response_channel = jch;
        % skip if all trials are bad
        if out.elecs(ich).all_bad(jch) == 1
            out.rejection_details(1).reject.keep(stim_channel,response_channel) = 0;
            out.rejection_details(2).reject.keep(stim_channel,response_channel) = 0;
            out.rejection_details(1).reject.all_bad(stim_channel,response_channel) = 1;
            out.rejection_details(2).reject.all_bad(stim_channel,response_channel) = 1;
            continue
        end
   
        % Get the baseline
        baseline = mean(eeg(1:stim_idx-idx_before_stim));
      
        % Get the eeg in the stim time
        stim_eeg = abs(eeg(temp_stim_idx(1):temp_stim_idx(2))-baseline);
        non_abs_stim_eeg = eeg(temp_stim_idx(1):temp_stim_idx(2))-baseline;
        
        % Get the eeg in the n1 and n2 time
        n1_eeg = eeg(temp_n1_idx(1):temp_n1_idx(2));
        n2_eeg = eeg(temp_n2_idx(1):temp_n2_idx(2));
        both_eeg = eeg(temp_n1_idx(1):temp_n2_idx(2));
        loose_n1_eeg = eeg(temp_loose_n1_idx(1):temp_loose_n1_idx(2));
        
        % subtract baseline
        n1_eeg_abs = abs(n1_eeg-baseline);
        n2_eeg_abs = abs(n2_eeg-baseline);
        both_eeg_abs = abs(both_eeg);
        
        % Get sd of baseline
        baseline_sd = std(eeg(1:stim_idx-idx_before_stim));

        % convert n1_eeg_abs to z score
        n1_z_score = n1_eeg_abs/baseline_sd;
        n2_z_score = n2_eeg_abs/baseline_sd;
        both_z_score = both_eeg_abs/baseline_sd;
        %}
        
        % OPTION A: FIRST PEAK 
        %% find the identity of the peaks
        choose_qualifier_MinLat = 0;
        n1_type_MinDist = 3; % DNE to start 
        [pks1_MinDist,locs1_MinDist] = findpeaks(n1_z_score,'MinPeakDistance',5e-3*fs);
        [peak1_MinDist,I1_MinDist]= min(locs1_MinDist);
        n1_peak_idx_MinDist = round(locs1_MinDist(I1_MinDist));
        n1_peak_MinDist = pks1_MinDist(I1_MinDist);
        I1_MinDist;
        if I1_MinDist>=0.01
            n1_peak_MinDist = n1_peak_MinDist;
            n1_peak_idx_MinDist = n1_peak_idx_MinDist;
        else
            n1_peak_MinDist = nan;
            n1_peak_idx_MinDist = nan;
        end
        if isempty(n1_peak_MinDist)
            n1_peak_MinDist = nan;
            n1_peak_idx_MinDist = nan;
        end
         
        if ~isnan(n1_peak_idx_MinDist)
            %n = 888
            n1_peak_idx_MinDist = n1_peak_idx_MinDist + temp_n1_idx(1) - 1 - stim_idx - 1;
            n1_objective_idx_MinDist = n1_peak_idx_MinDist + stim_start;
            n1_objective_idx_MinDist;
            n1_y_MinDist = out.elecs(ich).detrend_filt_avgs(n1_objective_idx_MinDist,jch);
            % n1_y_MinDist
            if n1_y_MinDist ==abs(n1_y_MinDist)
                n1_type_MinDist = 1; %maximum
            end
            if n1_y_MinDist ~= abs(n1_y_MinDist)
                n1_type_MinDist = 0; %minimum
            end
        end
        %n1_type_MinDist
        choose_qualifier_MinLat = choose_qualifier_MinLat + n1_peak_MinDist;
         
        n2_control_var_MinDist = 0;
        n2_type_MinDist = 4; % DNE to start 
        [pks2_MinDist,locs2_MinDist] = findpeaks(n2_z_score,'MinPeakDistance',5e-3*fs);
    
        while n2_control_var_MinDist ==0
            [n2_peak_MinDist,I2a_MinDist] = max(pks2_MinDist); % find the biggest
            %n2_peak_MinDist
            n2_peak_idx_MinDist = round(locs2_MinDist(I2a_MinDist));
            % n2_peak_idx_MinDist
            %isempty(n2_peak_MinDist)
            if isempty(n2_peak_MinDist)
                n2_peak_MinDist = nan;
                n2_peak_idx_MinDist = nan;
                n2_control_var_MinDist = 1;
            end
            if n2_peak_MinDist <2
                n2_peak_MinDist = nan;
                n2_peak_idx_MinDist = nan;
                n2_control_var_MinDist = 1;
            end
            if isnan(n2_peak_MinDist)
                n2_control_var_MinDist = 1;
                n2_peak_MinDist = nan;
                n2_peak_idx_MinDist = nan;
            end
            if ~isnan(n2_peak_idx_MinDist)
                n2_peak_idx_MinDist = n2_peak_idx_MinDist + temp_n2_idx(1) - 1 - stim_idx - 1;
                n2_objective_idx_MinDist = n2_peak_idx_MinDist + stim_start;
                n2_y_MinDist = out.elecs(ich).detrend_filt_avgs(n2_objective_idx_MinDist,jch);
                if n2_y_MinDist == abs(n2_y_MinDist)
                    n2_type_MinDist = 1; %maximum
                end
                if n2_y_MinDist ~= abs(n2_y_MinDist)
                    n2_type_MinDist = 0; %minimum
                end
                % n1_type_MinDist
                % n2_type_MinDist
                if n1_type_MinDist == n2_type_MinDist 
                    n2_control_var_MinDist = 1;
                end
                if n1_type_MinDist ~= n2_type_MinDist
                    pks2_MinDist(I2a_MinDist) = nan;
                end
            end
        end
        choose_qualifier_MinLat = choose_qualifier_MinLat + n2_peak_MinDist;

        % OPTION B: MAX PEAK 
        %% find the identity of the peaks
        choose_qualifier_MaxAmp = 0;
        n1_type_MaxAmp = 3; % DNE to start
        [pks1_MaxAmp,locs1_MaxAmp] = findpeaks(n1_z_score,'MinPeakDistance',5e-3*fs);
        [n1_peak_MaxAmp,I1_MaxAmp] = max(pks1_MaxAmp); % find the biggest
        n1_peak_idx_MaxAmp = round(locs1_MaxAmp(I1_MaxAmp));
        if I1_MaxAmp>=0.01
            n1_peak_MaxAmp = n1_peak_MaxAmp;
            n1_peak_idx_MaxAmp = n1_peak_idx_MaxAmp;
        else
            n1_peak_MaxAmp = nan;
            n1_peak_idx_MaxAmp = nan;
        end
        if isempty(n1_peak_MaxAmp)
            n1_peak_MaxAmp = nan;
            n1_peak_idx_MaxAmp = nan;
        end
        
        if ~isnan(n1_peak_idx_MaxAmp)
            n1_peak_idx_MaxAmp = n1_peak_idx_MaxAmp + temp_n1_idx(1) - 1 - stim_idx - 1;
            n1_objective_idx_MaxAmp = n1_peak_idx_MaxAmp + stim_start;
            n1_y_MaxAmp = out.elecs(ich).detrend_filt_avgs(n1_objective_idx_MaxAmp,jch);

            if n1_y_MaxAmp ==abs(n1_y_MaxAmp)
                n1_type_MaxAmp = 1; %maximum
            end
            if n1_y_MaxAmp ~= abs(n1_y_MaxAmp)
                n1_type_MaxAmp = 0; %minimum
            end
        end
        
        choose_qualifier_MaxAmp = choose_qualifier_MaxAmp + n1_peak_MaxAmp;

        n2_control_var_MaxAmp = 0;
        [pks2_MaxAmp,locs2_MaxAmp] = findpeaks(n2_z_score,'MinPeakDistance',5e-3*fs);
        n2_type_MaxAmp = 4; % DNE to start

        while n2_control_var_MaxAmp ==0
            [n2_peak_MaxAmp,I2a_MaxAmp] = max(pks2_MaxAmp); % find the biggest
            n2_peak_idx_MaxAmp = round(locs2_MaxAmp(I2a_MaxAmp));
            if isempty(n2_peak_MaxAmp)
                n2_peak_MaxAmp = nan;
                n2_peak_idx_MaxAmp = nan;
                n2_control_var_MaxAmp = 1;
            end
            if n2_peak_MaxAmp <2
                n2_peak_MaxAmp = nan;
                n2_peak_idx_MaxAmp = nan;
                n2_control_var_MaxAmp = 1;
            end
            if isnan(n2_peak_MaxAmp)
                n2_control_var_MaxAmp = 1;
                n2_peak_MaxAmp = nan;
                n2_peak_idx_MaxAmp = nan;
            end
            if ~isnan(n2_peak_idx_MaxAmp)
                n2_peak_idx_MaxAmp = n2_peak_idx_MaxAmp + temp_n2_idx(1) - 1 - stim_idx - 1;
                n2_objective_idx_MaxAmp = n2_peak_idx_MaxAmp + stim_start;
                n2_y_MaxAmp = out.elecs(ich).detrend_filt_avgs(n2_objective_idx_MaxAmp,jch);
                if n2_y_MaxAmp == abs(n2_y_MaxAmp)
                    n2_type_MaxAmp = 1; %maximum
                end
                if n2_y_MaxAmp ~= abs(n2_y_MaxAmp)
                    n2_type_MaxAmp = 0; %minimum
                end
                if n1_type_MaxAmp == n2_type_MaxAmp 
                    n2_control_var_MaxAmp = 1;
                end
                if n1_type_MaxAmp ~= n2_type_MaxAmp
                    pks2_MaxAmp(I2a_MaxAmp) = nan;
                end
                % if abs(n2_y_MaxAmp) < 0.5*(abs(n1_y_MaxAmp))
                %     n2_peak_MaxAmp = nan;
                %     n2_peak_idx_MaxAmp = nan;
                % end
            end
        end
        choose_qualifier_MaxAmp = choose_qualifier_MaxAmp + n2_peak_MaxAmp;

        if choose_qualifier_MinLat >= choose_qualifier_MaxAmp || n2_peak_idx_MaxAmp - n1_peak_idx_MaxAmp < 45
            n1_peak = n1_peak_MinDist;
            n1_peak_idx = n1_peak_idx_MinDist;
            n2_peak = n2_peak_MinDist;
            n2_peak_idx = n2_peak_idx_MinDist;
        else
            n1_peak = n1_peak_MaxAmp;
            n1_peak_idx = n1_peak_idx_MaxAmp;
            n2_peak = n2_peak_MaxAmp;
            n2_peak_idx = n2_peak_idx_MaxAmp;
        end
        
        % n1_peak
        % n1_peak_idx
        % n2_peak
        % n2_peak_idx

         % redefine idx relative to time after stim
        eeg_rel_peak_idx = n1_peak_idx + temp_n1_idx(1) - 1;
        %n1_peak_idx = n1_peak_idx + temp_n1_idx(1) - 1 - stim_idx - 1;
        %n2_peak_idx = n2_peak_idx + temp_n2_idx(1) - 1 - stim_idx - 1;
        
        if ~isnan(n1_peak_idx) % && ~isnan(n2_peak_idx)   
            % store   
            n1(jch,:) = [n1_peak,n1_peak_idx,nan,baseline_sd];
            n2(jch,:) = [n2_peak,n2_peak_idx,nan,baseline_sd];
        end
        % if isnan(n1_peak_idx) || isnan(n2_peak_idx)
        %     n1(jch,:) = [nan,nan];
        %     n2(jch,:) = [nan,nan];
        %     out.rejection_details(1).reject.keep(stim_channel,response_channel) = 0;
        %     out.rejection_details(2).reject.keep(stim_channel,response_channel) = 0;
        %     out.rejection_details(1).reject.no_both(stim_channel,response_channel)= 1;
        %     out.rejection_details(2).reject.no_both(stim_channel,response_channel) = 1;
        % end
        
        % if 0
        %     figure
        %     plot(eeg)
        %     hold on
        %     plot([temp_n1_idx(1) temp_n1_idx(1)],ylim)
        %     plot([temp_n1_idx(2) temp_n1_idx(2)],ylim)
        %     plot(xlim,[baseline baseline])
        % end
        
        %% Do various things to reject likely artifact

        % % 1:
        % % If sum of abs value in stim period is above a certain threshold
        % % relative to sum of abs value in n1 period, throw out n1
        % if sum(stim_eeg) > rel_thresh * sum(n1_eeg_abs)
        %     n1(jch,:) = [nan nan nan nan];
        %     out.rejection_details(1).reject.keep(stim_channel,response_channel) = 0;
        %     out.rejection_details(2).reject.keep(stim_channel,response_channel) = 0;
        % end
        % 
        % % 2:
        % % If anything too big in whole period, throw it out
        % if max(abs(eeg(temp_stim_idx(1):temp_n2_idx(end))-nanmedian(eeg))) > 1e3
        %     n1(jch,:) = [nan nan nan nan];
        %     n2(jch,:) = [nan nan nan nan]; 
        %     out.rejection_details(1).reject.keep(stim_channel,response_channel) = 0;
        %     out.rejection_details(2).reject.keep(stim_channel,response_channel) = 0;
        % end
        % 
        % 
        % % 3:
        % % If the EEG signal in the N1 period crosses a line connecting its
        % % first and last point more than twice, throw it out
        % n_crossings = count_crossings(loose_n1_eeg,baseline);
        % 
        % if n_crossings > max_crossings
        %     n1(jch,:) = [nan nan nan nan];
        %     n2(jch,:) = [nan nan nan nan];
        %     out.rejection_details(1).reject.keep(stim_channel,response_channel) = 0;
        %     out.rejection_details(2).reject.keep(stim_channel,response_channel) = 0;
        % end
        % %}
        % 
        % % 4:
        % % If no return to "baseline" between stim and N1, throw it out
        % 
        % return_to_baseline_before = 0;
        % signed_tight_stim_eeg = eeg(temp_tight_stim(1):temp_tight_stim(2))-baseline;
        % 
        % if ~isnan(n1_peak_idx)
        % 
        %     % if N1 above baseline
        %     if eeg(eeg_rel_peak_idx) - baseline > 0
        %         % Then look at max stim
        %         [max_stim,stim_max_idx] = max(signed_tight_stim_eeg);
        %         stim_height = max_stim - baseline;
        %         n1_height = eeg(eeg_rel_peak_idx) - baseline;
        %     else
        %         % Look at min stim
        %         [max_stim,stim_max_idx] = min(signed_tight_stim_eeg);
        %         stim_height = baseline - max_stim;
        %         n1_height =  baseline - eeg(eeg_rel_peak_idx);
        %     end
        %     stim_max_idx = stim_max_idx + temp_tight_stim(1) - 1;
        % 
        %     % Only invoke this rule if the height of the stim artifact
        %     % larger than height of n1
        %     if  stim_height > n1_height
        % 
        %         % If there's no part in between close to baseline
        %         bl_range = [baseline-1*baseline_sd,baseline+1*baseline_sd];
        % 
        %         if eeg(eeg_rel_peak_idx) - baseline > 0
        %              % check if it gets below the upper baseline range
        %             if any(eeg(stim_max_idx:eeg_rel_peak_idx) < bl_range(2))
        %                 return_to_baseline_before = 1;
        %             end
        %         else
        %             % check if it gets above the lower baseline range
        %             if any(eeg(stim_max_idx:eeg_rel_peak_idx) > bl_range(1))
        %                 return_to_baseline_before = 1;
        %             end
        %         end
        % 
        % 
        %         if 0
        %             figure
        %             plot(eeg)
        %             hold on
        %             plot(stim_max_idx,eeg(stim_max_idx),'o')
        %             plot(eeg_rel_peak_idx,eeg(eeg_rel_peak_idx),'o')
        %             plot(xlim,[bl_range(1) bl_range(1)])
        %             plot(xlim,[bl_range(2) bl_range(2)])
        %             if return_to_baseline_before
        %                 title('Ok')
        %             else
        %                 title('Not ok')
        %             end
        %         end
        % 
        %         if ~return_to_baseline_before
        %             n1(jch,:) = [nan nan nan nan];
        %             n2(jch,:) = [nan nan nan nan]; 
        %             out.rejection_details(1).reject.keep(stim_channel,response_channel) = 0;
        %             out.rejection_details(2).reject.keep(stim_channel,response_channel) = 0;
        %         end
        % 
        %     end
        % 
        % end
        % %}
        % 
        % % 5:
        % % if no return to baseline after N1 peak in a certain amount of
        % % time, throw it out
        % if ~isnan(n1_peak_idx)
        %     time_to_return_to_bl = 100e-3; % 50 ms
        %     idx_to_return_to_bl = eeg_rel_peak_idx+round(time_to_return_to_bl * fs);
        %     bl_range = [baseline-1.5*baseline_sd,baseline+1.5*baseline_sd];
        %     returns_to_baseline_after = 0;
        % 
        %     % if N1 above baseline
        %     if eeg(eeg_rel_peak_idx) - baseline > 0
        % 
        %         % check if it gets below the upper baseline range
        %         if any(eeg(eeg_rel_peak_idx:idx_to_return_to_bl) < bl_range(2))
        %             returns_to_baseline_after = 1;
        %         end
        % 
        %     else
        % 
        %         % check if it gets above the lower baseline range
        %         if any(eeg(eeg_rel_peak_idx:idx_to_return_to_bl) > bl_range(1))
        %             returns_to_baseline_after = 1;
        %         end
        % 
        %     end
        % 
        %     if ~returns_to_baseline_after
        %         n1(jch,:) = [nan nan nan nan];
        %         n2(jch,:) = [nan nan nan nan];
        %         out.rejection_details(1).reject.keep(stim_channel,response_channel) = 0;
        %         out.rejection_details(2).reject.keep(stim_channel,response_channel) = 0;
        %     end
        % 
        %     if 0
        %         figure
        %         plot(eeg)
        %         hold on
        %         plot([eeg_rel_peak_idx eeg_rel_peak_idx],ylim)
        %         plot([idx_to_return_to_bl...
        %             idx_to_return_to_bl],ylim)
        %         plot(xlim,[bl_range(1) bl_range(1)])
        %         plot(xlim,[bl_range(2) bl_range(2)])
        %         if returns_to_baseline_after
        %             title('Ok');
        %         else
        %             title('Not ok');
        %         end
        %     end
        % 
        % end

    end
    
    % add to struct
    out.elecs(ich).N1 = n1;
    out.elecs(ich).N2 = n2;

end
%% 
NumRows = size(out.elecs,2);
for ich = 1:NumRows
    if isempty(out.elecs(ich).arts)
        ppp = 9;
        out.elecs(ich).N1 = [];
        out.elecs(ich).N2 = [];
    end
end
%% 

function time = convert_idx_to_time(idx,times)

    time = linspace(times(1),times(2),length(idx));

end
end



