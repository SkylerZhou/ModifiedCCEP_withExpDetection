%% calcaulte spearman correlation between amplitude and latency. 

%% loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
patientNewOut_dir = locations.patientNewOut_dir;
patientOriout_dir = locations.patientOriOut_dir;
ptT = readtable([data_folder,'master_pt_list.xlsx']);

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

% to store mean and std of correlations across all electrodes for each patient 
n1_mean = zeros(1,num_patient);
n1_std = zeros(1,num_patient);
n2_mean = zeros(1,num_patient);
n2_std = zeros(1,num_patient);
all_n1_corrs = []; % To store all N1 correlations
all_n2_corrs = []; % To store all N2 correlations
idx_patient = 1;
%

for n = 1:num_patient

    which_version = 'new_pipeline';
    
    % obtain patient file 
    if strcmp(which_version, 'new_pipeline')
        patient_file = fullfile(patientNewOut_dir, patient_files(n));
        temp = load(patient_file);
        out = temp.new_out;
    else
        patient_file = fullfile(patientOriout_dir, patient_files(n));
        temp = load(patient_file);
        out = temp.pt_out;
    end

    % adjust amp and lat for n1&n2 to NaN if they are rejected in
    % rejection_details
    out = SZ_adjust_network_to_remove_rejects(out); 
    
    %% loop over elecs
    % to store correlations between the amplitude and latency of each electrode for
    % each patient 
    num_elecs = size(out.elecs, 2);
    non_empty_count = sum(~cellfun(@isempty, {out.elecs.arts})); % number of channel that started a stimulation 
    n1_corr_vector = zeros(1, non_empty_count);
    n2_corr_vector = zeros(1, non_empty_count);
    idx = 1;
    %

    %% for each stimulating electrode, what is the spearman correlations between its responding electrodes' amplitude and latencies 
    % --> should yield one correlation for each stimulating electrode 
    for ich = 1:num_elecs
        if isempty(out.elecs(ich).arts), continue; end

        % Calculate Spearman's correlation between N1 amplitude and latency
        n1_amp = out.elecs(ich).n1_adj(:,1);
        n1_lat = out.elecs(ich).n1_adj(:,2);
        n1_corr = corr(n1_amp, n1_lat,'rows','pairwise','type','spearman'); % NaN values will be ignored during pariwise 
        n1_corr_vector(idx) = n1_corr;
        
        % Calculate Spearman's correlation between N2 amplitude and latency
        n2_amp = out.elecs(ich).n2_adj(:,1);
        n2_lat = out.elecs(ich).n2_adj(:,2);
        n2_corr = corr(n2_amp, n2_lat,'rows','pairwise','type','spearman');
        n2_corr_vector(idx) = n2_corr;
    
        idx = idx + 1;
    end
    
    %% for each patient with num_elecs of stimulating electrode, what is the mean and std  
    n1_mean(idx_patient) = mean(n1_corr_vector,'omitnan');
    n1_std(idx_patient) = std(n1_corr_vector,'omitnan');
    n2_mean(idx_patient) = mean(n2_corr_vector','omitnan');
    n2_std(idx_patient) = std(n2_corr_vector,'omitnan');
    all_n1_corrs = [all_n1_corrs, n1_corr_vector];
    all_n2_corrs = [all_n2_corrs, n2_corr_vector];

    idx_patient = idx_patient + 1;
end

%% if the correlation is significantly negative
[h_n1, p_n1, ci_n1, stats_n1] = ttest(all_n1_corrs, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N1 correlations: ', num2str(p_n1)]);

[h_n2, p_n2, ci_n2, stats_n2] = ttest(all_n2_corrs, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N2 correlations: ', num2str(p_n2)]);



%% forest plot
figure;
tiledlayout(1,2);

% Number of patients
num_patient = length(n1_mean);

%% Plot for N1
nexttile;
hold on;

for i = 1:num_patient
    % Calculate the lower and upper bounds
    lower_bound = n1_mean(i) - n1_std(i);
    upper_bound = n1_mean(i) + n1_std(i);
    
    % Plot horizontal line representing the standard deviation interval
    line([lower_bound, upper_bound], [i, i], 'Color', 'black', 'LineWidth', 1);
    
    % Plot a marker at the mean
    plot(n1_mean(i), i, 'bo', 'MarkerFaceColor', [0 0.4470 0.7410]);
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
    % Calculate the lower and upper bounds
    lower_bound = n2_mean(i) - n2_std(i);
    upper_bound = n2_mean(i) + n2_std(i);
    
    % Plot horizontal line representing the standard deviation interval
    line([lower_bound, upper_bound], [i, i], 'Color', 'black', 'LineWidth', 1);
    
    % Plot a marker at the mean
    plot(n2_mean(i), i, 'bo', 'MarkerFaceColor', [0 0.4470 0.7410]);
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





%--------------------------------------------------------------------------
%% correlation of all electrodes at patient level. 
% obtain patient out file 
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

all_n1_corrs = nan(1, num_patient);
all_n2_corrs = nan(1, num_patient);


for n = 1:num_patient
    
    which_version = 'new_pipeline';
    
    % obtain patient file 
    patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', which_version, patient_files(n));
    temp = load(patient_file);
    if strcmp(which_version, 'new_pipeline')
        out = temp.new_out;
    else
        out = temp.pt_out;
    end

    % adjust amp and lat to remove those that are rejected for n1&n2 
    out = SZ_adjust_network_to_remove_rejects(out); 

    % to store amp and lat across all electrodes for each patient 
    num_elecs = size(out.elecs, 2);
    n1_ampLat_all = []; % To store N1 amp and lat pairs for all the stimulating electrode
    n2_ampLat_all = []; % To store N2 amp and lat pairs ..
    

    for ich = 1:num_elecs
        if isempty(out.elecs(ich).arts), continue; end
        n1_ampLat = out.elecs(ich).n1_adj;
        n2_ampLat = out.elecs(ich).n2_adj;
        
        % store the processed n1_ampLat 
        n1_ampLat = transpose(n1_ampLat);
        n1_ampLat_all = [n1_ampLat_all, n1_ampLat];

        % store the processed n2_ampLat 
        n2_ampLat = transpose(n2_ampLat);
        n2_ampLat_all = [n2_ampLat_all, n2_ampLat];
    end


    % calculate correlation between amplitude and latency all N1 and N2 
    n1_amp = n1_ampLat_all(1,:)';
    n1_lat = n1_ampLat_all(2,:)';
    n1_corr = corr(n1_amp, n1_lat,'rows','pairwise','type','spearman'); % NaN values will be ignored during pariwise 

    n2_amp = n2_ampLat_all(1,:)';
    n2_lat = n2_ampLat_all(2,:)';
    n2_corr = corr(n2_amp, n2_lat,'rows','pairwise','type','spearman'); % NaN values will be ignored during pariwise 

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
xlabel('Correlation beteween Amplitude and Latency');
title(sprintf('N1 Spearman Correlations Coefficient (p-value: %.4e)', p_n1));

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
xlabel('Correlation between Amplitdue and Latency');
title(sprintf('N2 Spearman Correlations Coefficient (p-value: %.4e)', p_n2));

grid on;
hold off;
%--------------------------------------------------------------------------





%--------------------------------------------------------------------------
%% scatter plot of amp and lat in patient level. 
% obtain patient out file 
num_patient = 20;
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(num_patient));


for n = num_patient:num_patient
    
    which_version = 'new_pipeline';
    
    % obtain patient file 
    patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', which_version, patient_files(n));
    temp = load(patient_file);
    if strcmp(which_version, 'new_pipeline')
        out = temp.new_out;
    else
        out = temp.pt_out;
    end

    % adjust amp and lat for n1&n2 
    out = SZ_adjust_network_to_remove_rejects(out); 

    % to store amp and lat across all electrodes for each patient 
    %n1_ampLat_oneElecs = []; % To store N1 amp and lat pairs for each stimulating electrode
    %n2_ampLat_oneElecs = []; % To store N2 amp and lat pairs ..
    num_elecs = size(out.elecs, 2);
    n1_ampLat_all = []; % To store N1 amp and lat pairs for all the stimulating electrode
    n2_ampLat_all = []; % To store N2 amp and lat pairs ..

    
    for ich = 1:num_elecs
        if isempty(out.elecs(ich).arts), continue; end
        n1_ampLat = out.elecs(ich).n1_adj;
        n2_ampLat = out.elecs(ich).n2_adj;
        
        % store the processed n1_ampLat 
        n1_ampLat = transpose(n1_ampLat);
        n1_ampLat_all = [n1_ampLat_all, n1_ampLat];
     
        % store the processed n2_ampLat 
        n2_ampLat = transpose(n2_ampLat);
        n2_ampLat_all = [n2_ampLat_all, n2_ampLat];
    end
end

%% scatter plot
scatter(n1_ampLat_all(2,:), n1_ampLat_all(1,:),...
    'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Set axes properties
xlim([0, 50]);
ylim([0, 150]);
xlabel('Latency (ms)');
ylabel('Amplitude');
title(sprintf('Amplitude vs. Latency for Patient %s', patient_ids));
