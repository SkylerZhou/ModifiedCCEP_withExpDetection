%% covert to csv with columns as patient, stim_channel, response_channel, elecs_dist, n1_amp, n1_lat
function [T,csv_name] = mat2csv(out, csv_folder)

%% extract N1, stim, and response channels info 
num_elecs = size(out.elecs, 2);
n1_all = [];
stim_all = [];
resp_all = [];
    
for ich = 1:num_elecs
    % if the current elecs is not a stimulating elecs, skip it
    %if isempty(out.elecs(ich).arts), continue; end

    % if exist distance info, extract n1 amp and lat, and distance
    if isfield(out.other, 'elecs_dist') 
        % cp the distance info as the third column to n1_adj
        %out.elecs(ich).n1_adj(:,3) = out.other.elecs_dist(:,ich);
        % extract n1 amp and lat, and distance
        %n1 = out.elecs(ich).n1_adj;
        n1 = out.other.elecs_dist(:,ich);
        n1_all = [n1_all; n1];
    %else
        % else extract only n1 amp and lat
        %n1 = out.elecs(ich).n1_adj;
        %n1_all = [n1_all; n1];
    end

    % extract stim and response labels 
    stim = repmat(out.bipolar_labels(ich), size(out.bipolar_labels,1), 1);
    stim_all = [stim_all; stim];
    resp = out.bipolar_labels(1: size(out.bipolar_labels,1));
    resp_all = [resp_all; resp];
end


%% create cell arrays containing patient, stim_channel, response_channel, elecs_dist
num_rows = size(stim_all, 1);
num_cols = 4; % patient, stim_channel, response_channel, elecs_dist
new_out = cell(num_rows, num_cols);
    
% add patient id to col 1
patient_col = repmat({out.name}, num_rows, 1);
new_out(:,1) = patient_col;
    
% add stim_channel and resp_channel to col 2 and 3
new_out(:,2) = stim_all;
new_out(:,3) = resp_all;

% append elecs_dist to col 4 
n1_all = num2cell(n1_all);
if size(n1_all, 2) == 1
    % If distance information is available
    %new_out(:,5) = n1_all(:,1); % n1 amp
    %new_out(:,6) = n1_all(:,2); % n1 lat
    new_out(:,4) = n1_all; % elecs_dist
else
    % If distance information is not available
    %new_out(:,5) = n1_all(:,1); % n1 amp
    %new_out(:,6) = n1_all(:,2); % n1 lat
    new_out(:,4) = num2cell(NaN(num_rows,1)); % Placeholder for elecs_dist
end



%% save the new output cell array file as csv
varNames = {'PatientID', 'StimChannel', 'RespChannel', 'ElecsDist'};
T = cell2table(new_out, 'VariableNames', varNames);
csv_name = fullfile(csv_folder, [out.name, '.csv']);


  