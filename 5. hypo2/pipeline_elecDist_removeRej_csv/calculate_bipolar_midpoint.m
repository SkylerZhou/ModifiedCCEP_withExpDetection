%% this function calculate elucidean distance betwwen electrodes and store the distance matrix as elecs_distance in out.other?
function out = calculate_bipolar_midpoint(out)

% Mount to remote server and define file path
locations = cceps_files;
mount_dir = locations.mount_dir;
patient_id = out.name;
coor_filename = string(strcat(patient_id, '.csv'));


%% Attempt to load the coordinates file
try
    coor_file = readtable(string(fullfile(mount_dir, coor_filename)));
    coor_file.name = string(coor_file.labels);  % Convert label column to 'name' for consistency
catch ME
    fprintf('Error finding the patient electrode coordinates file %s: %s\n', coor_filename, ME.message);
    out.other.bipolar_coor = table();  % Return empty table if failed
    return;
end


%% Convert bipolar labels to string for later split with -
bi_elecs = cell2table(out.bipolar_labels, 'VariableNames', {'name'});
bi_elecs.names = cellfun(@(x) string(x), bi_elecs.name, 'UniformOutput', false);
catch_idx = cellfun(@(x) isempty(x) || ~ischar(x), bi_elecs.name);
bi_elecs.names(catch_idx) = {""};
bi_elecs = removevars(bi_elecs, 'name');



%% Create a map from monopolar electrode name to XYZ coordinates
names_to_xyz = containers.Map(coor_file.name, num2cell([coor_file.mm_x, coor_file.mm_y, coor_file.mm_z], 2));

% Initialize output holders (preallocate to preserve order)
n = height(bi_elecs);
bipolar_labels = bi_elecs.names;  % same order
bipolar_coords = NaN(n, 3);      % default to NaN

% Loop through bipolar electrodes
for i = 1:n
    label = string(bi_elecs.names(i));  % e.g., "LA1-LA2"

    if label == ""
        continue;  % leave as NaN
    end

    parts = strsplit(label, '-');
    elec1 = parts{1};
    elec2 = parts{2};

    if isKey(names_to_xyz, elec1) && isKey(names_to_xyz, elec2)
        coord1 = names_to_xyz(elec1);
        coord2 = names_to_xyz(elec2);
        midpoint = (coord1 + coord2) / 2;
        bipolar_coords(i, :) = midpoint;  % insert at correct index
    end
end


%% Store output in a table
out.other.bipolar_coor = table( ...
    bipolar_labels, ...
    bipolar_coords(:,1), bipolar_coords(:,2), bipolar_coords(:,3), ...
    'VariableNames', {'label', 'mm_x', 'mm_y', 'mm_z'});

