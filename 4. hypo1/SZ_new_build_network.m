% to see gui, need to click on a specific x&y, press enter, click outside
% the figure, and finally press spacebar

function out = SZ_new_build_network(out,do_gui,save_name)

if ~exist('do_gui','var'), do_gui = 0; end
if ~exist('save_name','var'), save_name = ''; end

%% Parameters
thresh_amp = 4.5;
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

% Loop over n1 and n2
for w = 1:length(wavs)

    which = wavs{w};
    % change the number of rows for A
    A = nan(nchs,nchs);

    %% initialize rejection details

    % ---- edited on Mar 19, 2025 by skyler -----------
    %{
    details.thresh = thresh_amp;
    details.which = which;
    details.reject.sig_avg = nan(length(elecs),length(elecs));
    details.reject.pre_thresh = nan(length(elecs),length(elecs));
    details.reject.at_thresh = nan(length(elecs),length(elecs));
    details.reject.keep = nan(length(elecs),length(elecs));
    details.reject.no_n1 = nan(length(elecs),length(elecs));
    details.reject.no_both = nan(length(elecs),length(elecs));
    details.reject.exp = nan(length(elecs),length(elecs));
    details.reject.empty = nan(length(elecs),length(elecs));
    %}
    % -----------------------------------------------

    for ich = 1:length(elecs)

        if isempty(elecs(ich).arts), continue; end
        
        % ---- edited on Mar 19, 2025 by skyler ------------ 
        % Add peak amplitudes to the array
        % - which should use n1_adj instead? - or use n1 to see their amp as
        % well. 
        % - only the ich that correspond to the sz onset can be selected?
        arr = elecs(ich).(which);
        A(ich,:) = arr(:,1); % get the amplitude 
        % --------------------------------------------------


        % ---- edited on Mar 19, 2025 by skyler ------------ 
        %{
        all_bad = logical(elecs(ich).all_bad);
        details.reject.sig_avg(ich,:) = all_bad;
        details.reject.pre_thresh(ich,:) = isnan(elecs(ich).(which)(:,1)) & ~all_bad;
        details.reject.at_thresh(ich,:) = elecs(ich).(which)(:,1) < thresh_amp;
        details.reject.keep(ich,:) = elecs(ich).(which)(:,1) >= thresh_amp;
        all_nans = (sum(~isnan(elecs(ich).avg),1) == 0)';
        details.reject.sig_avg(ich,:) = all_nans;
        details.reject.pre_thresh(ich,:) = isnan(elecs(ich).(which)(:,1)) & ~all_nans;
        details.reject.at_thresh(ich,:) = elecs(ich).(which)(:,1) < thresh_amp;
        details.reject.keep(ich,:) = elecs(ich).(which)(:,1) >= thresh_amp;
        %}
        % ---------------------------------------------
    end
    
    % ---- edited on Mar 19, 2025 by skyler ------------ 
    % Add details to array
    % out.rejection_details(w) = details;
    % ---------------------------------------------

    %% Remove ignore chs
    stim_chs = nansum(A,2) > 0;
    response_chs = keep_chs;
    A(:,~response_chs) = nan;
    A = A'; % transverse!!
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
    out.network(w).A0 = A0;

end


if do_gui == 1
%% PLot
figure
set(gcf,'position',[1 11 1400 900])
show_network_no_fig(out,1,0,0,save_name)
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
    drawnow;
    pause;
    close(gcf);
end
end
%}

end