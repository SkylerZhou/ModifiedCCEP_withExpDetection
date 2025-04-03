function ccep_assessment(out, num_subplots, pt_id, sz_id, stim_resp)


%% Parameters for Plotting 
n1_time = [11e-3 50e-3];     
n2_time = [50e-3 300e-3];     
zoom_times = [-300e-3 300e-3];  
zoom_factor = 2;            
which_n = 1;                 
which = out.rejection_details(which_n).which;

sig_overlap_stims = stim_resp.ccep_sig_overlap_stims;
sig_overlap_resps = stim_resp.ccep_sig_overlap_resps;
sig_only_stims = stim_resp.ccep_sig_only_stims;
sig_only_resps = stim_resp.ccep_sig_only_resps;
nonsig_overlap_stims = stim_resp.ccep_nonsig_overlap_stims;
nonsig_overlap_resps = stim_resp.ccep_nonsig_overlap_resps;

%% Create a New Figure
fig = figure('Color', 'w', 'Visible', 'on');

% Determine Grid Layout for Subplots
n_per_line = ceil(sqrt(num_subplots));   % Number of columns in the grid
n_lines = ceil(num_subplots / n_per_line); % Number of rows in the grid
t = tiledlayout(n_lines, n_per_line, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Determine the Number of Pairs in Each Group
num_sig_overlap = numel(sig_overlap_stims); 
num_sig_only = numel(sig_only_stims);    
num_nonsig_overlap = numel(nonsig_overlap_stims);    

% Initialize counters for each group
sig_overlap_idx = 1;
sig_only_idx    = 1;
nonsig_overlap_idx = 1;

%% Loop Through Each Subplot
for i = 1:num_subplots        
    ax = nexttile(t);
    ax.FontSize = 6;        
    hold(ax, 'on');   

    % Decide whether this subplot corresponds to an overlapping or an only pair.
    if sig_overlap_idx <= num_sig_overlap
          
        % --- Plot for a Sig CCEP Overlap Pair ---    
        % get the row and col idx for out 
        stim_ch = sig_overlap_stims{sig_overlap_idx};
        resp_ch = sig_overlap_resps{sig_overlap_idx};
        row = find(strcmp(out.bipolar_labels, stim_ch));
        col = find(strcmp(out.bipolar_labels, resp_ch));
            
        % extract the wave data from out 
        avg = out.elecs(row).detrend_filt_avgs(:,col);
        times = out.elecs(row).times;
        eeg_times = convert_indices_to_times(1:length(avg),out.other.stim.fs,times(1));
        wav =  out.elecs(row).(which)(col,:);
                    
        stim_idx = out.elecs(row).stim_idx;
        wav_idx = wav(2)+stim_idx+1;
        wav_time = convert_indices_to_times(wav_idx,out.other.stim.fs,times(1));
        n1_idx = floor(n1_time*out.other.stim.fs);
        n2_idx = floor(n2_time*out.other.stim.fs);
        temp_n1_idx = n1_idx + stim_idx - 1;
        temp_n2_idx = n2_idx + stim_idx - 1;

        % plot with a red line for overlapping pairs
        plot(ax, eeg_times, avg, 'Color', [0.8500 0.3250 0.0980], 'linewidth', 1.5);
        title_str = sprintf('Stim %s, Resp %s', num2str(stim_ch), num2str(resp_ch));
        hold on;

        % annotate N1 location and z-score 
        if (out.elecs(row).N1(col,1))~=0
            x = (out.elecs(row).N1(col,2)/out.other.stim.fs);
            x_indx = round(out.elecs(row).N1(col,2)+stim_idx+1); 
            y = out.elecs(row).detrend_filt_avgs(x_indx,col);
            plot(x,y,'bX','markersize',8,'linewidth',1.5);
            text(wav_time+0.01,avg(round(wav_idx)),sprintf('%s z-score: %1.1f', which,wav(1)), 'fontsize',6)
        end
        hold on
                    
        % annotate N2 location               
        if ~isnan(out.elecs(row).N2(col,2))               
            x = (out.elecs(row).N2(col,2)/out.other.stim.fs);                               
            x_indx = out.elecs(row).N2(col,2)+stim_idx+1;                
            y = out.elecs(row).detrend_filt_avgs(x_indx,col);                
            plot(x,y,'rX','markersize',8,'linewidth',1.5);                
        end

        % set xlim and ylim
        xlim([zoom_times(1) zoom_times(2)]);
        height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2))-median(avg)));
        if ~any(isnan(avg))
            ylim([median(avg)-zoom_factor*height,median(avg)+zoom_factor*height]);
        end  
        plot([0 0],ylim,'k--');
            
        % next plot 
        sig_overlap_idx = sig_overlap_idx + 1;


    elseif sig_overlap_idx > num_sig_overlap && sig_only_idx <= num_sig_only

        % --- Plot for a Sig CCEP Only Pair ---
        % get the row and col idx for out 
        stim_ch = sig_only_stims{sig_only_idx};
        resp_ch = sig_only_resps{sig_only_idx};
        row = find(strcmp(out.bipolar_labels, stim_ch));
        col = find(strcmp(out.bipolar_labels, resp_ch));

        % extract the wave data from out 
        avg = out.elecs(row).detrend_filt_avgs(:,col);
        times = out.elecs(row).times;
        eeg_times = convert_indices_to_times(1:length(avg),out.other.stim.fs,times(1));
        wav =  out.elecs(row).(which)(col,:);
                    
        stim_idx = out.elecs(row).stim_idx;
        wav_idx = wav(2)+stim_idx+1;
        wav_time = convert_indices_to_times(wav_idx,out.other.stim.fs,times(1));
        n1_idx = floor(n1_time*out.other.stim.fs);
        n2_idx = floor(n2_time*out.other.stim.fs);
        temp_n1_idx = n1_idx + stim_idx - 1;
        temp_n2_idx = n2_idx + stim_idx - 1;

        % plot with a blue line for only pairs
        plot(ax, eeg_times, avg, 'Color', [0.4660 0.6740 0.1880], 'linewidth', 1.5);
        title_str = sprintf('Stim %s, Resp %s', num2str(stim_ch), num2str(resp_ch));
        hold on;

        % annotate N1 location and z-score 
        if (out.elecs(row).N1(col,1))~=0
            x = (out.elecs(row).N1(col,2)/out.other.stim.fs);         
            x_indx = round(out.elecs(row).N1(col,2)+stim_idx+1); 
            y = out.elecs(row).detrend_filt_avgs(x_indx,col);
            plot(x,y,'bX','markersize',8,'linewidth',1.5);
            text(wav_time+0.01,avg(round(wav_idx)),sprintf('%s z-score: %1.1f', which,wav(1)), 'fontsize',6)
        end
        hold on
                
        % annotate N2 location               
        if ~isnan(out.elecs(row).N2(col,2))
            x = (out.elecs(row).N2(col,2)/out.other.stim.fs);                
            x_indx = out.elecs(row).N2(col,2)+stim_idx+1;
            y = out.elecs(row).detrend_filt_avgs(x_indx,col);
            plot(x,y,'rX','markersize',8,'linewidth',1.5);                
        end

        % set xlim and ylim
        xlim([zoom_times(1) zoom_times(2)]);
        height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2))-median(avg)));
        if ~any(isnan(avg))
            ylim([median(avg)-zoom_factor*height,median(avg)+zoom_factor*height]);
        end  
        plot([0 0],ylim,'k--');
            
        % next plot 
        sig_only_idx = sig_only_idx + 1;

    else

        % --- Plot for a NonSig CCEP Overlap Pair ---
        % get the row and col idx for out 
        stim_ch = nonsig_overlap_stims{nonsig_overlap_idx};
        resp_ch = nonsig_overlap_resps{nonsig_overlap_idx};
        row = find(strcmp(out.bipolar_labels, stim_ch));
        col = find(strcmp(out.bipolar_labels, resp_ch));

        % extract the wave data from out 
        avg = out.elecs(row).detrend_filt_avgs(:,col);
        times = out.elecs(row).times;
        eeg_times = convert_indices_to_times(1:length(avg),out.other.stim.fs,times(1));
        wav =  out.elecs(row).(which)(col,:);
                    
        stim_idx = out.elecs(row).stim_idx;
        wav_idx = wav(2)+stim_idx+1;
        wav_time = convert_indices_to_times(wav_idx,out.other.stim.fs,times(1));
        n1_idx = floor(n1_time*out.other.stim.fs);
        n2_idx = floor(n2_time*out.other.stim.fs);
        temp_n1_idx = n1_idx + stim_idx - 1;
        temp_n2_idx = n2_idx + stim_idx - 1;

        % plot with a blue line for only pairs
        plot(ax, eeg_times, avg, 'Color', [0 0.4470 0.7410], 'linewidth', 1.5);
        title_str = sprintf('Stim %s, Resp %s', num2str(stim_ch), num2str(resp_ch));
        hold on;

        % annotate N1 location and z-score 
        if (out.elecs(row).N1(col,1))~=0
            x = (out.elecs(row).N1(col,2)/out.other.stim.fs);         
            x_indx = round(out.elecs(row).N1(col,2)+stim_idx+1); 
            y = out.elecs(row).detrend_filt_avgs(x_indx,col);
            plot(x,y,'bX','markersize',8,'linewidth',1.5);
            text(wav_time+0.01,avg(round(wav_idx)),sprintf('%s z-score: %1.1f', which,wav(1)), 'fontsize',6)
        end
        hold on
                
        % annotate N2 location               
        if ~isnan(out.elecs(row).N2(col,2))
            x = (out.elecs(row).N2(col,2)/out.other.stim.fs);                
            x_indx = out.elecs(row).N2(col,2)+stim_idx+1;
            y = out.elecs(row).detrend_filt_avgs(x_indx,col);
            plot(x,y,'rX','markersize',8,'linewidth',1.5);                
        end

        % set xlim and ylim
        xlim([zoom_times(1) zoom_times(2)]);
        height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2))-median(avg)));
        if ~any(isnan(avg))
            ylim([median(avg)-zoom_factor*height,median(avg)+zoom_factor*height]);
        end  
        plot([0 0],ylim,'k--');
            
        % next plot 
        nonsig_overlap_idx = nonsig_overlap_idx + 1;

    end
        
    % Add title and labels for each subplot
    title(ax, title_str, 'Interpreter', 'none', 'FontSize', 6);
    grid(ax, 'on');
end

%% Remove Unused Subplots 
totalTiles = n_lines * n_per_line;
if totalTiles > num_subplots
    for j = (num_subplots+1):totalTiles
        axEmpty = nexttile(t);
        axis(axEmpty, 'off');
    end
end

% Add a Main Figure Title
sgtitle_str = sprintf('CCEP Assessment for Patient %s with Seizure ID %s.', pt_id, sz_id);
sgtitle(fig, sgtitle_str, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');

% Adjust Figure Layout
set(fig, 'Units', 'normalized', 'Position', [0.05 0.05 1.0 1.0]);

% Save the Figure as PNG to a Specific Folder
output_folder = '/Users/zhouzican/Documents/CNT/DSOSD/SZ_hypo1_data/ccep_assessment'; 
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
filename = fullfile(output_folder, sprintf('CCEP_Assessment_Patient_%s_SZ_%s.png', pt_id, sz_id));
saveas(fig, filename);

end
