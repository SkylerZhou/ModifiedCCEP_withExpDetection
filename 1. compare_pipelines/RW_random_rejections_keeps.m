function RW_random_rejections_keeps(out)

%% Parameters
pretty = 0;
n_to_plot = 25; % how many total to show
n_per_line = 5;
n_lines = 5;
n1_time = [11e-3 50e-3];
n2_time = [50e-3 300e-3];
zoom_times = [-300e-3 300e-3];
zoom_factor = 2;
which_n = 1;

%% Get various path locations
locations = cceps_files; % Need to make a file pointing to you own path
pwfile = locations.pwfile;
loginname = locations.loginname;
script_folder = locations.script_folder;
results_folder = locations.results_folder;
thirdOut_dir = locations.thirdOut_dir;

% add paths
addpath(genpath(script_folder));
if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end
name = out.name;
out_folder = [thirdOut_dir,'validation/',name,'/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% Pick intracranial chs with bipolar signal
keep_chs = get_chs_to_ignore(out.bipolar_labels); % get_chs_to_ignore removes bad chs

%% Get rejection details arrays
thresh = out.rejection_details(which_n).thresh;
which = out.rejection_details(which_n).which;

sig_avg = out.rejection_details(which_n).reject.sig_avg;
pre_thresh = out.rejection_details(which_n).reject.pre_thresh;
at_thresh = out.rejection_details(which_n).reject.at_thresh;
% sz commented out 2025/04/06
%no_both = out.rejection_details(which_n).reject.no_both; 
keep = out.rejection_details(which_n).reject.keep;
exp = out.rejection_details(which_n).reject.exp;
ignore_ch = out.rejection_details(which_n).reject.ignore_ch;

% sz commented out 2025/04/06
% any_reject = sig_avg == 1| pre_thresh == 1 | at_thresh == 1 | no_both == 1 | exp ==1 ;
any_reject = sig_avg == 1| pre_thresh == 1 | at_thresh == 1 | exp ==1 | ignore_ch==1 ;

% Calculate total numbers
nkeep = sum(keep(:) == 1);
nreject = sum(any_reject(:) == 1);
nunstim = sum(isnan(keep(:)));

% if nunstim+nreject+nkeep ~= size(keep,1)*size(keep,1)
%      error('numbers do not add up');
% end

% Loop through rejection types
for j = 1:2
    if j == 1
        thing = keep;
        cat = 'New Keep';
    else
        thing = any_reject;
        cat = 'Reject Any';
    end
    
    meet_criteria = find(thing==1);

    % sz comment out on 04/14/2025 with the addition of rejection_details.ignore_ch
    % Restrict to those on keep chs
    %[row,col] = ind2sub(size(keep),meet_criteria);  % obtain the linear indices of all the keeps 
    %meet_criteria(keep_chs(row) == false) = []; % filter out the keeps that are recorded with the bad channels(not sure)
    %col(keep_chs(row) == false) = [];
    %meet_criteria(keep_chs(col) == false) = [];
    
    % Initialize figure
    figure
    set(gcf,'position',[100 100 1200 1000])
    t = tiledlayout(n_lines,n_per_line,'padding','compact','tilespacing','compact');
    
 
    % Pick a random N
    to_plot = randsample(meet_criteria,n_to_plot); % randomly select 25 from the coordinates of all the keeps/rejects that meets the criteria 
    
    % Loop through these
    for i=1:size(to_plot)
        
        ind = to_plot(i);
        
        % convert this to row and column 
        [row,col] = ind2sub(size(keep),ind);  % obtain the coordinate of a certain linear indices of a randomly selected keep/reject, etc 7299 --> [101,60]
        
        row;
        col;
        
        % get why it was rejected
        why = nan;
        if j == 2 
            if sig_avg(row,col) == 1
                why = 'averaging';
            end
            if pre_thresh(row,col) == 1
                % if ~isnan(why)
                %     error('what');
                % end
                why = 'artifact';
            end
            if at_thresh(row,col) == 1
                % if ~isnan(why)
                %     error('what');
                % end
                why = 'threshold';
            end
            % sz commented out 2025/04/06
            %if no_both(row,col) == 1
                %if isnan(why)
                    %why = 'no both';
                %end
            %end
            if exp(row,col) == 1
                if isnan(why)
                    why = 'exponential';
                end
            end
            if ignore_ch(row,col) == 1
                if isnan(why)
                    why = 'ignore_ch';
                end
            end           
        end
        
        % Get the waveform
        try
            avg = out.elecs(row).detrend_filt_avgs(:,col);
        catch ME
            % If there's an error, print the row, col, and size of the array
            fprintf('Error at row: %d, col: %d\n', row, col);
            fprintf('Size of out.elecs(%d).detrend_filt_avgs: %s\n', row, mat2str(size(out.elecs(row).detrend_filt_avgs)));
            fprintf('nkeep: %d, nreject: %d\n', nkeep, nreject);
            rethrow(ME); % Rethrow the error to stop execution
        end
            
        %avg = out.elecs(row).detrend_filt_avgs(:,col);
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
        
        % plot the keeps 
        if j==1
            if ~isnan(out.elecs(row).N1(col,1))
                % Plot the detrend_filt_avgs
                nexttile
                plot(eeg_times,avg,'k','linewidth',2);
                hold on
                
                % annotate N1 location and z-score 
                if (out.elecs(row).N1(col,1))~=0
                    x = (out.elecs(row).N1(col,2)/out.other.stim.fs);
                    x_indx = round(out.elecs(row).N1(col,2)+stim_idx+1);
                    y = out.elecs(row).detrend_filt_avgs(x_indx,col);
                    plot(x,y,'bX','markersize',15,'linewidth',4);
                    if ~pretty
                    text(wav_time+0.01,avg(round(wav_idx)),sprintf('%s z-score: %1.1f',...
                        which,wav(1)), 'fontsize',9)
                    end
                end
                hold on
                
                % annotate N2 location
                if ~isnan(out.elecs(row).N2(col,2))
                    x = (out.elecs(row).N2(col,2)/out.other.stim.fs);
                    x_indx = out.elecs(row).N2(col,2)+stim_idx+1;
                    y = out.elecs(row).detrend_filt_avgs(x_indx,col);
                    plot(x,y,'rX','markersize',15,'linewidth',4);
                end
        
                %xlim([eeg_times(1) eeg_times(end)])
                xlim([zoom_times(1) zoom_times(2)]);
                
                % Zoom in (in the y-dimension) around the maximal point in the N1
                % time period
                height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2))-median(avg)));
                if ~any(isnan(avg))
                    ylim([median(avg)-zoom_factor*height,median(avg)+zoom_factor*height]);
                end
                
                % annoatate stim and response electrodes, as well as stim
                % start time 
                labels = out.bipolar_labels;
                stim_label = labels{row};
                resp_label = labels{col};  
                stim_start = out.other.periods(row).start_time;
                pause(0.1)
                xl = xlim;
                yl = ylim;
                if ~pretty
                    text(xl(1),yl(2),sprintf('Stim: %s\nResponse: %s\nStart: %.2f s',stim_label,resp_label,stim_start),...
                        'horizontalalignment','left',...
                        'verticalalignment','top','fontsize',9);
                end
                plot([0 0],ylim,'k--');
                set(gca,'fontsize',9)
                if pretty
                    yticklabels([])
                    %xticklabels([])
                    xtl = xticklabels;
                    xtlc = cellfun(@(x) sprintf('%s s',x),xtl,'uniformoutput',false);
                    %xlabel('Time (s)')
                    xticklabels(xtlc)
                end
            end
        end



        % plot the rejects
        if j == 2
            %if ~isnan(out.elecs(row).N1(col,1))
                % Plot the detrend_filt_avg
                nexttile
                plot(eeg_times,avg,'k','linewidth',2);
                hold on
                
                % annoatate N1 location and z-score
                if (out.elecs(row).N1(col,1))~=0 && ~isnan(out.elecs(row).N1(col,1))
                    x = (out.elecs(row).N1(col,2)/out.other.stim.fs);
                    x_indx = round(out.elecs(row).N1(col,2)+stim_idx+1);
                    y = out.elecs(row).detrend_filt_avgs(x_indx,col);
                    plot(x,y,'bX','markersize',15,'linewidth',4);
                    if ~pretty
                    text(wav_time+0.01,avg(round(wav_idx)),sprintf('%s z-score: %1.1f',...
                        which,wav(1)), 'fontsize',9)
                    end
                end
                hold on

                % annotate N2 location
                if ~isnan(out.elecs(row).N2(col,2))
                    x = (out.elecs(row).N2(col,2)/out.other.stim.fs);
                    x_indx = out.elecs(row).N2(col,2)+stim_idx+1
                    y = out.elecs(row).detrend_filt_avgs(x_indx,col);
                    plot(x,y,'rX','markersize',15,'linewidth',4);
                end
        
                % set xlim and ylim
                %xlim([eeg_times(1) eeg_times(end)])
                xlim([zoom_times(1) zoom_times(2)]);
                
                % Zoom in (in the y-dimension) around the maximal point in the N1
                % time period
                height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2))-median(avg)));
                if ~any(isnan(avg))
                    ylim([median(avg)-zoom_factor*height,median(avg)+zoom_factor*height]);
                end
                

                % annotate stim and response electroes, as well as
                % stim_start time 
                labels = out.bipolar_labels;
                stim_label = labels{row};
                resp_label = labels{col};
                stim_start = out.other.periods(row).start_time;
                pause(0.1)
                xl = xlim;
                yl = ylim;
                if ~pretty
                    text(xl(1),yl(2),sprintf('Stim: %s\nResponse: %s\nStart: %.2f s',stim_label,resp_label,stim_start),...
                        'horizontalalignment','left',...
                        'verticalalignment','top','fontsize',9);
                end
                plot([0 0],ylim,'k--');
                set(gca,'fontsize',9)
                if pretty
                    yticklabels([])
                    %xticklabels([])
                    xtl = xticklabels;
                    xtlc = cellfun(@(x) sprintf('%s s',x),xtl,'uniformoutput',false);
                    %xlabel('Time (s)')
                    xticklabels(xtlc)
                end
                title(why)
            %end
        end
    end
    

    % save the figures 
    if pretty == 0
        title(t,sprintf('%s %s z-score threshold %1.1f',cat,which,thresh));
    end
    
    if pretty
        fname = sprintf('%s_%sthresh_%d_pretty.png',cat,which,thresh);
    else
        fname = sprintf('%s_%sthresh_%d.png',cat,which,thresh);
    end
    
    print(gcf,[out_folder,fname],'-dpng');
    
end