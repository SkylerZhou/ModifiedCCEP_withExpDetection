%% this function calculate elucidean distance betwwen electrodes and store the distance matrix as elecs_distance in out.other?
function out = add_elecs_distance(out)


% mount to remote server
locations = cceps_files;
mount_dir = locations.mount_dir;
patient_id = out.name;
coordinates_filename = string(strcat(patient_id, '.csv'));


% read the metadata file containing the electrode from the mounted dir
% names and xyz coordinates 
% load coordinates_file if it is avilable in the remote server
try
    coordinates_file = readtable(string(strcat(mount_dir, coordinates_filename)));
catch ME 
    % if error occurs, suggests the coordinates info for the
    % corresponding HUP ID does not exist 
    fprintf('Error finding the patient electrode coordinates file %s: %s\n', coordinates_filename, ME.message);
end


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
