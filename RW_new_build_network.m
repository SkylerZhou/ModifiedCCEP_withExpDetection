% Construct a connectivity network from the CCEP using the filtered kept
% signals.
% Stimulation given at 512 Hz. Recording/Sampling rate is 2 * 512 Hz = 1024 Hz
% Total time duration is -0.5 + 0.8 = 1.3 secs
% Number of samples recorded during this time window = 1.3 sec * 1024 Hz =
% 1332, which can be seen in size(out.elecs(ich).avg and detrend_filt_avgs)
function out = RW_new_build_network(out,do_gui) 

if ~exist('do_gui','var'), do_gui = 0; end

%% Parameters
thresh_amp = 6.0;
thresh_exp = 0.6;
peak_end_time = 0.3;
wavs = {'N1','N2'};

%% Get various path locations
locations = cceps_files; % Need to make a file pointing to you own path
script_folder = locations.script_folder;
results_folder = locations.results_folder;

% add paths
addpath(genpath(script_folder));
if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

%% Basic info
elecs = out.elecs;
chLabels = out.chLabels;
nchs = length(chLabels);
keep_chs = get_chs_to_ignore(chLabels);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% skyler's code: add in goodness of fit in out.elecs

% coef range for the fitted exponential function
aLower = -10;
aUpper = 10;
bLower = -30;
bUpper = 0;

% Pre-allocate the gof matrix with nan
for row = 1:nchs
    if size(out.elecs(row).avg, 1) >= 1
        out.elecs(row).gof = nan(1, nchs); 
    end
end

for row=1:nchs
    times = out.elecs(row).times;
    stim_idx = out.elecs(row).stim_idx;
    for col=1:nchs
        if size(out.elecs(row).avg,1) >= 1
            if out.elecs(row).N1(col,1) > 0
            
                % if stim-response pair, extract the timeseries which starts after peak_idx and end before 0.3 sec.
                peak_start_index = out.elecs(row).N1(col,2) + stim_idx;  % the index of 0 sec is 512
                peak_end_index = stim_idx + floor(out.other.stim.fs * peak_end_time);  
                
                eeg_times = convert_indices_to_times(peak_start_index:peak_end_index,out.other.stim.fs,times(1));
                after_peak_avg = out.elecs(row).detrend_filt_avgs(peak_start_index:peak_end_index,col); 
                
                % standarize the after_peak_avg between [-1,1]
                mean_val = mean(after_peak_avg);
                std_val = std(after_peak_avg);
                if std_val ~= 0
                    after_peak_avg_std = (after_peak_avg - mean_val) / std_val;
                else
                    after_peak_avg_std = zeros(size(after_peak_avg)); % Handle constant signal
                end

                % Define cutoff frequency and filter parameters
                cutoff_freq = 50;  
                Fs = 2 * out.other.stim.fs;
                [b, a] = butter(4, cutoff_freq/(Fs/2), 'low');  % Fs is the sampling frequency
                after_peak_avg_std = filtfilt(b, a, after_peak_avg_std);

                % examine goodness of fit 
                fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'Lower', [aLower, bLower], ...
                        'Upper', [aUpper, bUpper]);
                [f, gof] = fit(eeg_times',after_peak_avg_std(:),'exp1', fitOptions);
                out.elecs(row).gof(:,col) = gof.adjrsquare;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Loop over n1 and n2. analyzes the two types of signal separately
for w = 1:length(wavs)

    which = wavs{w};
    A = nan(nchs,nchs);  % initialize a matrix A with NaN. 

    %% initialize rejection details
    details.thresh = thresh_amp;
    details.which = which;
    % initialize matrices within out.rejection_details.reject to track the
    % rejection criteria
    details.reject.sig_avg = nan(length(elecs),length(elecs)); % track signal that are identified as bad based on averaging criteria
    details.reject.pre_thresh = nan(length(elecs),length(elecs)); % identifies signals below the threshold before applying thresh_amp?
    details.reject.at_thresh = nan(length(elecs),length(elecs)); % identifies if signals meet the thresh_amp
    details.reject.keep = nan(length(elecs),length(elecs));
    details.reject.no_n1 = nan(length(elecs),length(elecs));
    details.reject.no_both = nan(length(elecs),length(elecs));
    %details.reject.deriv = nan(length(elecs),length(elecs));
    details.reject.exp = nan(length(elecs),length(elecs));

    % for each of N1 and N2 wave, loop over each electrode
    for ich = 1:length(elecs)
        
        % if no stimulation from the electrode, skip it 
        if isempty(elecs(ich).arts), continue; end
        ich;
        if size(out.elecs(ich).arts,1)~=0
            arr = elecs(ich).(which); % elecs(ich).('N1') stores the peak amplitude and latencies of all the other electrodes after stimulation at this ich electrode 
        
            % Add peak amplitudes (which is arr(:,1)) to the row corresponding to the stimulating electrode in A matrix
            A(ich,:) = arr(:,1);
        
            all_bad = logical(elecs(ich).all_bad);
            details.reject.sig_avg(ich,:) = all_bad;
            details.reject.pre_thresh(ich,:) = isnan(elecs(ich).(which)(:,1)) & ~all_bad;
            details.reject.at_thresh(ich,:) = elecs(ich).(which)(:,1) < thresh_amp;
            details.reject.keep(ich,:) = elecs(ich).(which)(:,1) >= thresh_amp;
            details.reject.no_n1(ich,:) = isnan(elecs(ich).(which)(:,1));
            
            %{
            all_nans = (sum(~isnan(elecs(ich).avg),1) == 0)';
            details.reject.sig_avg(ich,:) = all_nans;
            details.reject.pre_thresh(ich,:) = isnan(elecs(ich).(which)(:,1)) & ~all_nans;
            details.reject.at_thresh(ich,:) = elecs(ich).(which)(:,1) < thresh_amp;
            details.reject.keep(ich,:) = elecs(ich).(which)(:,1) >= thresh_amp;
            %}
        end
    end

    % Add details to array
    out.rejection_details(w) = details;


    %% Remove ignore chs
    stim_chs = nansum(A,2) > 0;
    response_chs = keep_chs;
    A(:,~response_chs) = nan;
    A = A';
    A0 = A;

    %% Normalize
    A(A0<thresh_amp) = 0;
    
    %% Add this to array
    if w == 1
        out.stim_chs = stim_chs;
        out.response_chs = response_chs;
    end
    
    out.network(w).which = which;
    out.network(w).A = A;

end

n = size(out.chLabels,1);
for ich=1:n
    for jch=1:n
        if size(out.elecs(ich).N1,1) >=2
            if ~isnan(out.elecs(ich).N1(jch,1)) && isnan(out.elecs(ich).N2(jch,1))   
                out.rejection_details(1).reject.keep(ich,jch) = 0;
                out.rejection_details(2).reject.keep(ich,jch) = 0;
                out.rejection_details(1).reject.no_both(ich,jch)= 1;
                out.rejection_details(2).reject.no_both(ich,jch) = 1;
            end
        end
        
        %{
        % set 1/10 of the max of the de-trended averages as the threshold
        % of slope 
        if size(out.elecs(ich).deriv,1) >= 1
            deriv_upper = max(out.elecs(ich).deriv(:,jch));
            deriv_upper_limit = (1/10)*(max(out.elecs(ich).detrend_filt_avgs(:,jch)));
            if deriv_upper >=deriv_upper_limit
                out.rejection_details(1).reject.keep(ich,jch) = 0;
                out.rejection_details(2).reject.keep(ich,jch) = 0;
                out.rejection_details(1).reject.deriv(ich,jch) = 1;
                out.rejection_details(2).reject.deriv(ich,jch) = 1;
            end
        end
        %}
 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % skyler's code: evaluate goodness of fit, reject cceps with high
        % gof to exponential curve
        if size(out.elecs(ich).avg,1) >= 1
            if out.elecs(ich).N1(jch,1) > 0
                if out.elecs(ich).gof(:,jch) >= thresh_exp
                    out.rejection_details(1).reject.keep(ich,jch) = 0;
                    out.rejection_details(2).reject.keep(ich,jch) = 0;
                    out.rejection_details(1).reject.exp(ich,jch) = 1;
                    out.rejection_details(2).reject.exp(ich,jch) = 1;
                end
            end           
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

    end
end




%% Convert electrode labels to anatomic locations
%{
response_labels = chLabels(response_chs);
stim_labels = chLabels(stim_chs);
mean_positions_response = 1:length(response_labels);
mean_positions_stim = 1:length(stim_labels);
edge_positions_response = [];
edge_positions_stim = []; 
    

ch_info.response_chs = response_chs;
ch_info.stim_chs = stim_chs;
ch_info.stim_pos = mean_positions_stim;
ch_info.response_pos = mean_positions_response;
ch_info.stim_labels = stim_labels;
ch_info.response_labels = response_labels;
ch_info.normalize = normalize;
ch_info.response_edges = edge_positions_response;
ch_info.stim_edges = edge_positions_stim;
ch_info.waveform = which;
%}


%{
in_degree = nansum(A,2);
[in_degree,I] = sort(in_degree,'descend');
in_degree_chs = chs(I);

out_degree = nansum(A,1);
[out_degree,I] = sort(out_degree,'descend');
out_degree_chs = chs(I);

if isempty(ana)
    all_labels = chLabels;
else
    all_labels = ana(chs);
end
ana_word = justify_labels(all_labels,'none');

if normalize == 1 || normalize == 0
fprintf('\nThe highest in-degree channels (note normalization!) are:\n');
for i = 1:min(10,length(in_degree))
    fprintf('%s (%s) (in-degree = %1.1f)\n',...
        chLabels{in_degree_chs(i)},ana_word{in_degree_chs(i)},in_degree(i));
end
end

if normalize == 2 || normalize == 0
fprintf('\nThe highest out-degree channels (note normalization!) are:\n');
for i = 1:min(10,length(out_degree))
    fprintf('%s (%s) (out-degree = %1.1f)\n',...
        chLabels{out_degree_chs(i)},ana_word{out_degree_chs(i)},out_degree(i));
end
end
%}


% n = 246;
% for ich=1:n
%     for jch=1:n
%         if size(out.elecs(ich).N1,1)~=0
%             if isnan(out.elecs(ich).N1(jch,1))
%                 if out.rejection_details(1).reject.keep(ich,jch) ~=1 
%                    if out.rejection_details(1).reject.at_thresh(ich,jch)~=1
%                         if out.rejection_details(1).reject.pre_thresh(ich,jch)~=1
%                             if out.rejection_details(1).reject.sig_avg(ich,jch)~=1
%                                 out.rejection_details(1).reject.no_n1(ich,jch)=1;
%                                 out.rejection_details(2).reject.no_n1(ich,jch)=1;
%                             end
%                         end
%                    end
%                 end
%             end
%         end
%     end
% end



if do_gui == 1
%% PLot
figure
set(gcf,'position',[1 11 1400 900])
show_network_no_fig(out,1,1,0)
stim_ch_idx = find(stim_chs);
response_ch_idx = find(response_chs);

while 1
    try
        [x,y] = ginput;
    catch
        break
    end
    if length(x) > 1, x = x(end); end
    if length(y) > 1, y = y(end); end
    figure
    set(gcf,'position',[215 385 1226 413])
    %tight_subplot(1,1,[0.01 0.01],[0.15 0.10],[.02 .02]);
    show_avg(out,stim_ch_idx(round(x)),response_ch_idx(round(y)),0,1)
    
    pause
    close(gcf)
end
end
%}

end