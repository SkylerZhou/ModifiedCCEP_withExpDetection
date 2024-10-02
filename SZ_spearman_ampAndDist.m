%% calculate spearman correlation between amplitude and electrode distance


%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
ptT = readtable([data_folder,'master_pt_list.xlsx']);

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

% array to store the correlations
all_n1_corrs = nan(1, num_patient);
all_n2_corrs = nan(1, num_patient);
%


%% loop over patients
for n = 1:num_patient
    

    % load patient out file 
    patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(n));
    temp = load(patient_file);
    out = temp.new_out;
    

    % try add_elecs_distance function; check if there are corresponding
    % coordinate files to build distance matrix for this patient 
    try 
        out = SZ_add_elecs_distance(out);
    catch ME 
        fprintf('%s\n', ME.message);
        continue
    end
    % adjust amp and lat to remove those that are rejected for n1&n2 
    out = SZ_adjust_network_to_remove_rejects(out); 


    % to store amp and lat across all electrodes for each patient 
    num_elecs = size(out.elecs, 2);
    n1_ampDist_all = []; % To store N1 amp and lat pairs for all the stimulating electrode
    n2_ampDist_all = []; % To store N2 amp and lat pairs ..
    

    for ich = 1:num_elecs
        if isempty(out.elecs(ich).arts), continue; end

        % cp the distance info as the thrid column to n1_adj&n2_adj
        out.elecs(ich).n1_adj(:,3) = out.other.elecs_dist(:,ich);
        out.elecs(ich).n2_adj(:,3) = out.other.elecs_dist(:,ich);

        n1_ampDist = out.elecs(ich).n1_adj(:,[1,3]);
        n2_ampDist = out.elecs(ich).n2_adj(:,[1,3]);
        
        % store the processed n1_ampLat 
        n1_ampDist = transpose(n1_ampDist);
        n1_ampDist_all = [n1_ampDist_all, n1_ampDist];

        % store the processed n2_ampLat 
        n2_ampDist = transpose(n2_ampDist);
        n2_ampDist_all = [n2_ampDist_all, n2_ampDist];
    end


    % calculate correlation between amplitude and latency all N1 and N2 
    n1_amp = n1_ampDist_all(1,:)';
    n1_dist = n1_ampDist_all(2,:)';
    n1_corr = corr(n1_amp, n1_dist,'rows','pairwise','type','spearman'); % NaN values will be ignored during pariwise 

    n2_amp = n2_ampDist_all(1,:)';
    n2_dist = n2_ampDist_all(2,:)';
    n2_corr = corr(n2_amp, n2_dist,'rows','pairwise','type','spearman'); % NaN values will be ignored during pariwise 

    all_n1_corrs(n) = n1_corr;
    all_n2_corrs(n) = n2_corr;
end


%% if the correlation is significantly negative
[h_n1, p_n1, ci_n1, stats_n1] = ttest(all_n1_corrs, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N1 correlations: ', num2str(p_n1)]);

[h_n2, p_n2, ci_n2, stats_n2] = ttest(all_n2_corrs, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N2 correlations: ', num2str(p_n2)]);


%% forest plot
figure;
tiledlayout(1,2);

%% Plot for N1
nexttile;
hold on;

for i = 1:num_patient   
    % Plot a marker at the mean
    plot(all_n1_corrs(i), i, 'bo', 'MarkerFaceColor', [0 0.4470 0.7410]);
end

% Draw vertical red line at x=0
xline(0, ':', 'LineWidth', 1);

% Set axes properties
xlim([-1, 1]);
ylim([0.5, num_patient + 0.5]);
set(gca, 'YTick', 1:num_patient, 'YTickLabel', patient_ids);
ylabel('Patient ID');
xlabel('Spearman Correlation');
title(sprintf('Forest Plot of N1 Correlations (p-value: %.4e)', p_n1));

grid on;
hold off;

%% Plot for N2
nexttile;
hold on;

for i = 1:num_patient
    % Plot a marker at the mean
    plot(all_n2_corrs(i), i, 'bo', 'MarkerFaceColor', [0 0.4470 0.7410]);
end

% Draw vertical red line at x=0
xline(0, ':', 'LineWidth', 1);

% Set axes properties
xlim([-1, 1]);
ylim([0.5, num_patient + 0.5]);
set(gca, 'YTick', 1:num_patient, 'YTickLabel', patient_ids);
ylabel('Patient ID');
xlabel('Spearman Correlation');
title(sprintf('Forest Plot of N2 Correlations (p-value: %.4e)', p_n2));

grid on;
hold off;
%--------------------------------------------------------------------------



