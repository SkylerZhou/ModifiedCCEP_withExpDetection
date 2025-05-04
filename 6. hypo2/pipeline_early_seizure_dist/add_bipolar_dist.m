%% this function calculate elucidean distance betwwen electrodes and store the distance matrix as elecs_distance in out.other?
function out = add_bipolar_dist(out)


%% mount to remote server
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
    out.other.bipolar_dist = table();  % Return empty table if failed
    return;
end


%% Load bipolar coordinate table
bipolar_table = out.other.bipolar_coor;

% Extract coordinates and check for valid entries
coords = [bipolar_table.mm_x, bipolar_table.mm_y, bipolar_table.mm_z];
valid_mask = all(~isnan(coords), 2);  % true for rows with valid XYZ

n = height(bipolar_table);
dist_matrix = NaN(n);  % initialize with NaN


%% Compute pairwise distances only for valid coordinate rows
for i = 1:n
    if ~valid_mask(i), continue; end
    for j = i:n
        if ~valid_mask(j), continue; end
        d = norm(coords(i,:) - coords(j,:));  % Euclidean distance
        dist_matrix(i,j) = d;
        dist_matrix(j,i) = d;  % symmetry
    end
end

% Store in output
out.other.bipolar_dist = array2table(dist_matrix);




%{
% delete on May 3, 2025 by sz
% check whether the elecs number and names in the coordinate_file and
% patient_file match; adjust electrodes and their corresponding xyz
% coordinates as nan if not matched
patient_elecs = cell2table(out.chLabels, 'VariableNames', {'name'});
patient_elecs.name = string(patient_elecs.name);
coordinates_file.name = string(coordinates_file.labels);

% have to re-order the output after outerjoin because outerjoin
% automatically sort by number 
[coor_file, rowsidx_left] = outerjoin(patient_elecs, coordinates_file, 'Keys', 'name', 'Type', 'left');
[~, sortinds] = sort(rowsidx_left);
coor_file = coor_file(sortinds,:);

% compute the pairwise euclidean distances between electrodes
coor_file = coor_file(:,3:5); % keep only xyz to compute distance
elecs_dist = pdist2(table2array(coor_file), table2array(coor_file), 'euclidean');

% create elecs_dist matrix under out.other
out.other.elecs_dist = elecs_dist;
%}
