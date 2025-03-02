%% calculate spearman correlation between amplitude and electrode distance on the electrode level


%% prep 
% load directories to loop over patients
locations = cceps_files;
data_folder = locations.data_folder;
patientNewOut_dir = locations.patientNewOut_dir;
ptT = readtable([data_folder,'master_pt_list.xlsx']);

num_patient = height(ptT);
patient_files = string(strcat(ptT.HUPID, '.mat'));
patient_ids = string(ptT.HUPID(1:num_patient));


n1_ampDist_corr_all = []; % To store the correlation of N1 amp and dist pairs for all electrodes
n2_ampDist_corr_all = []; % To store for N2
elec_counter = 0;
    

%% loop over patients
for n = 1:num_patient

    % load patient out file 
    patient_file = fullfile(patientNewOut_dir, patient_files(n));
    temp = load(patient_file);
    out = temp.out;
    
    % try add_elecs_distance function; check if there are corresponding
    % coordinate files to build distance matrix for this patient 
    try 
        out = SZ_add_elecs_distance(out);
    catch ME 
        fprintf('%s\n', ME.message);
        continue
    end
    % adjust amp and lat to remove those that are rejected for n1&n2 if not already did in SZ_runNew_compareOriNew.m
    % out = SZ_adjust_network_to_remove_rejects(out); 


    % to store amp and dist across all electrodes  
    num_elecs = size(out.elecs, 2);

    for ich = 1:num_elecs
        if isempty(out.elecs(ich).arts), continue; end

        temp_data = out.elecs(ich).n1_adj;
        temp_numPairs = sum(any(~isnan(temp_data), 2));
        elec_counter = elec_counter + temp_numPairs; 

        % cp the distance info as the third column to n1_adj&n2_adj
        out.elecs(ich).n1_adj(:,3) = out.other.elecs_dist(:,ich);
        out.elecs(ich).n2_adj(:,3) = out.other.elecs_dist(:,ich);

        % get n1_amp and n1_dist for each electrode; similarly obtain these
        % for n2
        n1_amp = out.elecs(ich).n1_adj(:,1);
        n1_dist = out.elecs(ich).n1_adj(:,3);
        n2_amp = out.elecs(ich).n2_adj(:,1);
        n2_dist = out.elecs(ich).n2_adj(:,3);

        % calculate the spearman's correlation between amp and dist for
        % this electrode
        n1_ampDist_corr = corr(n1_amp, n1_dist,'rows','pairwise','type','spearman');
        n2_ampDist_corr = corr(n2_amp, n2_dist,'rows','pairwise','type','spearman');

        % store the correlation of the single electrode
        n1_ampDist_corr_all = [n1_ampDist_corr_all, n1_ampDist_corr];
        n2_ampDist_corr_all = [n2_ampDist_corr_all, n2_ampDist_corr];      
    end
end

n1_ampDist_corr_avg = mean(n1_ampDist_corr_all,2, 'omitmissing');
n2_ampDist_corr_avg = mean(n2_ampDist_corr_all,2, 'omitmissing');
%





%--------------------------------------------------------------------------
%% scatter plot of amp and distance in patient level.

n1_ampDist_all = []; % To store N1 lat and dist pairs for all the stimulating electrode
n2_ampDist_all = []; % To store N2 lat and dist pairs ..

% loop over patients
for n = 1:num_patient

    % load patient out file 
    patient_file = fullfile(patientNewOut_dir, patient_files(n));
    temp = load(patient_file);
    out = temp.out;
    
    % try add_elecs_distance function; check if there are corresponding
    % coordinate files to build distance matrix for this patient 
    try 
        out = SZ_add_elecs_distance(out);
    catch ME 
        fprintf('%s\n', ME.message);
        continue
    end


    % to store amp and lat across all electrodes for each patient 
    num_elecs = size(out.elecs, 2);
    
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
end

numofelectpairs = sum(~isnan(n1_ampDist_all(1,:)));

%% scatter plot
% for n1
scatter(n1_ampDist_all(1,:), n1_ampDist_all(2,:),...
    'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Set axes properties
xlim([0, 220]);
ylim([0, 340]);
xlabel('Amplitude');
ylabel('Distance');
title('N1 Amplitude vs. Distance Scatter Plot');

txt = sprintf('rho = %f; #pts = 40; #elecs = %i', n1_ampDist_corr_avg, elec_counter);
text(100,150,txt)



% for n2
scatter(n2_ampDist_all(1,:), n2_ampDist_all(2,:),...
    'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Set axes properties
xlim([0, 200]);
ylim([0, 150]);
xlabel('Latency (ms)');
ylabel('Distance');
title(sprintf('Latency vs. Distance for Patient %s', patient_ids));




