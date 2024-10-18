%% calculate spearman correlation between amplitude and electrode distance


%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
patientNewOut_dir = locations.patientNewOut_dir;
ptT = readtable([data_folder,'master_pt_list.xlsx']);

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));

% array to store the correlations
all_n1_ampLat_corr = nan(1, num_patient);
all_n1_ampDist_corr = nan(1, num_patient);
all_n1_latDist_corr = nan(1, num_patient);

all_n2_ampLat_corr = nan(1, num_patient);
all_n2_ampDist_corr = nan(1, num_patient);
all_n2_latDist_corr = nan(1, num_patient);
%


%% loop over patients
for n = 1:1

    % load patient out file 
    patient_file = fullfile(patientNewOut_dir, patient_files(n));
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
    n1_ampLatDist_all = []; % To store N1 amp, lat, and dist pairs for all the stimulating electrode
    n2_ampLatDist_all = []; % To store N2 amp, lat, and dist pairs ..
    

    for ich = 1:num_elecs
        if isempty(out.elecs(ich).arts), continue; end

        % cp the distance info as the thrid column to n1_adj&n2_adj
        out.elecs(ich).n1_adj(:,3) = out.other.elecs_dist(:,ich);
        out.elecs(ich).n2_adj(:,3) = out.other.elecs_dist(:,ich);
        
        % store the processed n1_ampLatDist and n2_ampLatDist 
        n1_ampLatDist = transpose(out.elecs(ich).n1_adj);
        n2_ampLatDist = transpose(out.elecs(ich).n2_adj);

        n1_ampLatDist_all = [n1_ampLatDist_all, n1_ampLatDist];
        n2_ampLatDist_all = [n2_ampLatDist_all, n2_ampLatDist];        
    end


    % calculate correlation between amplitude and latency all N1 and N2 
    n1_amp = n1_ampLatDist_all(1,:)';
    n1_lat = n1_ampLatDist_all(2,:)';
    n1_dist = n1_ampLatDist_all(3,:)';
    n1_ampLat_corr = corr(n1_amp, n1_lat,'rows','pairwise','type','spearman'); % NaN values will be ignored during pariwise 
    n1_ampDist_corr = corr(n1_amp, n1_dist,'rows','pairwise','type','spearman');
    n1_latDist_corr = corr(n1_lat, n1_dist,'rows','pairwise','type','spearman');

    n2_amp = n2_ampLatDist_all(1,:)';
    n2_lat = n2_ampLatDist_all(2,:)';
    n2_dist = n2_ampLatDist_all(3,:)';
    n2_ampLat_corr = corr(n2_amp, n2_lat,'rows','pairwise','type','spearman'); % NaN values will be ignored during pairwise 
    n2_ampDist_corr = corr(n2_amp, n2_dist,'rows','pairwise','type','spearman');
    n2_latDist_corr = corr(n2_lat, n2_dist,'rows','pairwise','type','spearman');

    all_n1_ampLat_corr(n) = n1_ampLat_corr;
    all_n1_ampDist_corr(n) = n1_ampDist_corr;
    all_n1_latDist_corr(n) = n1_latDist_corr;

    all_n2_ampLat_corr(n) = n2_ampLat_corr;
    all_n2_ampDist_corr(n) = n2_ampDist_corr;
    all_n2_latDist_corr(n) = n2_latDist_corr;
end
%



%% if the correlation is significantly negative
[~, p_n1_ampLat, ~, ~] = ttest(all_n1_ampLat_corr, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N1 correlations between amplitude and latency: ', num2str(p_n1_ampLat)]);
[~, p_n1_ampDist, ~, ~] = ttest(all_n1_ampDist_corr, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N1 correlations between amplitude and distance: ', num2str(p_n1_ampDist)]);
[~, p_n1_latDist, ~, ~] = ttest(all_n1_latDist_corr, 0, 'Tail', 'right'); 
disp(['p-value for the one-sample t-test on N1 correlations between latency and distance: ', num2str(p_n1_latDist)]);

[~, p_n2_ampLat, ~, ~] = ttest(all_n2_ampLat_corr, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N2 correlations between amplitude and latency: ', num2str(p_n2_ampLat)]);
[~, p_n2_ampDist, ~, ~] = ttest(all_n2_ampDist_corr, 0, 'Tail', 'left'); 
disp(['p-value for the one-sample t-test on N2 correlations between amplitude and distance: ', num2str(p_n2_ampDist)]);
[~, p_n2_latDist, ~, ~] = ttest(all_n2_latDist_corr, 0, 'Tail', 'right'); 
disp(['p-value for the one-sample t-test on N2 correlations between latency and distance: ', num2str(p_n2_latDist)]);
%




%% forest plot
% has to manually change the following variables to plot for different
% graphs
%{
x_label = 'Amplitude and Latency';
toPlot_n1 = all_n1_ampLat_corr;
toPlot_n2 = all_n2_ampLat_corr;
p_n1 = p_n1_ampLat;
p_n2 = p_n2_ampLat;
%}

x_label = 'Amplitude and Distance';
toPlot_n1 = all_n1_ampDist_corr;
toPlot_n2 = all_n2_ampDist_corr;
p_n1 = p_n1_ampDist;
p_n2 = p_n2_ampDist;

%{
x_label = 'Latency and Distance';
toPlot_n1 = all_n1_latDist_corr;
toPlot_n2 = all_n2_latDist_corr;
p_n1 = p_n1_latDist;
p_n2 = p_n2_latDist;
%}
%


%% to plot
figure;
tiledlayout(1,2);
%


%% Plot for N1 correlations 
nexttile;
hold on;

for i = 1:num_patient   
    % Plot a marker at the mean
    plot(toPlot_n1(i), i, 'bo', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
end

% Draw vertical red line at x=0
xline(0, ':', 'LineWidth', 1);

% Set axes properties
xlim([-1, 1]);
ylim([0.5, num_patient + 0.5]);
set(gca, 'YTick', 1:num_patient, 'YTickLabel', patient_ids);
ylabel('Patient ID');
xlabel(sprintf('Correlation between %s', x_label));
title(sprintf('N1 Spearman Correlations Coefficienet (p-value: %.4e)', p_n1));

grid on;
hold off;
%


%% Plot for N2
nexttile;
hold on;

for i = 1:num_patient
    % Plot a marker at the mean
    plot(toPlot_n2(i), i, 'bo', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
end

% Draw vertical red line at x=0
xline(0, ':', 'LineWidth', 1);

% Set axes properties
xlim([-1, 1]);
ylim([0.5, num_patient + 0.5]);
set(gca, 'YTick', 1:num_patient, 'YTickLabel', patient_ids);
ylabel('Patient ID');
xlabel(sprintf('Correlation between %s', x_label));
title(sprintf('N2 Spearman Correlations Coefficient (p-value: %.4e)', p_n2));

grid on;
hold off;
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
%% scatter plot of lat and distance in patient level. 

%% prep
locations = cceps_files;
data_folder = locations.data_folder;
patientNewOut_dir = locations.patientNewOut_dir;
ptT = readtable([data_folder,'master_pt_list.xlsx']);

num_patient = 32;
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(num_patient));
%


%% obtain the patient file; get lat and distance 
for n = num_patient:num_patient
    
    which_version = 'new_pipeline';
    
    % obtain patient file 
    patient_file = fullfile(patientNewOut_dir, patient_files(n));
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
    n1_latDist_all = []; % To store N1 lat and dist pairs for all the stimulating electrode
    n2_latDist_all = []; % To store N2 lat and dist pairs ..

    
    for ich = 1:num_elecs
        if isempty(out.elecs(ich).arts), continue; end

        % cp the distance info as the thrid column to n1_adj&n2_adj
        out.elecs(ich).n1_adj(:,3) = out.other.elecs_dist(:,ich);
        out.elecs(ich).n2_adj(:,3) = out.other.elecs_dist(:,ich);

        n1_latDist = out.elecs(ich).n1_adj(:,[2,3]);
        n2_latDist = out.elecs(ich).n2_adj(:,[2,3]);
        
        % store the processed n1_ampLat 
        n1_latDist = transpose(n1_latDist);
        n1_latDist_all = [n1_latDist_all, n1_latDist];
     
        % store the processed n2_ampLat 
        n2_latDist = transpose(n2_latDist);
        n2_latDist_all = [n2_latDist_all, n2_latDist];
    end
end


%% scatter plot
% for n1
scatter(n1_latDist_all(1,:), n1_latDist_all(2,:),...
    'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Set axes properties
xlim([0, 50]);
ylim([0, 150]);
xlabel('Latency (ms)');
ylabel('Distance');
title(sprintf('Latency vs. Distance for Patient %s', patient_ids));


% for n2
scatter(n2_latDist_all(1,:), n2_latDist_all(2,:),...
    'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Set axes properties
xlim([0, 200]);
ylim([0, 150]);
xlabel('Latency (ms)');
ylabel('Distance');
title(sprintf('Latency vs. Distance for Patient %s', patient_ids));
