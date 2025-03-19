function show_network_no_fig(out,which,do_log,do_save,save_name)


%% Get various path locations
locations = cceps_files; % Need to make a file pointing to you own path
script_folder = locations.script_folder;
results_folder = locations.results_folder;
thirdOut_dir = locations.thirdOut_dir;
out_folder = [thirdOut_dir,'pretty_nets/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end


% setup
stim_chs = out.stim_chs;
response_chs = out.response_chs;
A = out.network(which).A;
C = out.name;
C = strsplit(C,'_');
name = C{1};


% remove channels not marked as stim/resp
A(~response_chs,:) = [];
A(:,~stim_chs) = [];

if do_log
    A = log(A);
    ltext = '(log scale)';
else
    ltext = '';
end

if which == 1
    wtext = 'N1';
elseif which == 2
    wtext = 'N2';
end


% plot the network 
imagesc(A);
hold on
xticklabels([])
yticklabels([])
xlabel('Stimulation site')
ylabel('Response site')
set(gca,'fontsize',20)
c = colorbar;
ylabel(c,[wtext,' z-score',ltext]);
title([name,' ',wtext])


% ------ edited on Mar 19, 2029 by skyler -------
stimidx = find(stim_chs);
respidx = find(response_chs);

% check the corresponding n1 amp from the stim channel's n1_adj.
for r = 1:length(respidx)
    for c = 1:length(stimidx)
        stim_idx = stimidx(c);
        resp_idx = respidx(r);
        % If the n1 amplitude for this stim-response pair is not nan,
        % overlay an asterisk on figure 
        if ~isnan(out.elecs(stim_idx).n1_adj(resp_idx,1))
            plot(c, r, 'm*', 'MarkerSize', 6, 'LineWidth', 1);
        end
    end
end
% ------------------------------------------------


if do_save
    print([out_folder,save_name,wtext],'-dpng')
    close(gcf)
end


end
