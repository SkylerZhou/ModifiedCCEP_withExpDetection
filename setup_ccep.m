% a global function to switch between different versions of CCEP
% run before using a specific version:
% setup_ccep('name_of_specific_version_used') 

function setup_ccep(version)
    % Define base directory where all versions and required folders are stored
    base_dir = '/Users/zhouzican/Documents/MATLAB/master_thesis/';
    
    % Define paths for different CCEP versions
    versions = struct( ...
        'CCEP', fullfile(base_dir, 'CCEP'), ...
        'CCEP_Detection', fullfile(base_dir, 'CCEP_Detection'), ...
        'CCEP_ExpDetection', fullfile(base_dir, 'CCEP_ExpDetection'), ...
        'CCEP_Revised', fullfile(base_dir, 'CCEP_Revised') ...
    );

    % Define paths that should always be retained
    permanent_paths = { ...
        fullfile(base_dir, 'ieeg_login'), ...
        fullfile(base_dir, 'ieeg-matlab-1.14.60') ...
    };

    % Remove all previously added CCEP version paths
    rmpath(genpath(base_dir));

    % Add back permanent paths
    for i = 1:length(permanent_paths)
        if isfolder(permanent_paths{i})
            addpath(genpath(permanent_paths{i}));
        else
            warning('Permanent path not found: %s', permanent_paths{i});
        end
    end
    
    % Add the selected CCEP version
    if isfield(versions, version)
        addpath(genpath(versions.(version)));
        disp(['Using CCEP version: ', version]);

        % Special case: If CCEP_ExpDetection is selected, also add CCEP_Revised but exclude cceps_files.m
        if strcmp(version, 'CCEP_ExpDetection')
            revised_path = versions.CCEP_Revised;

            % Get all subfolders of CCEP_Revised
            subfolders = strsplit(genpath(revised_path), pathsep);

            % Add each subfolder manually, but avoid adding the folder containing cceps_files.m
            for i = 1:length(subfolders)
                if exist(fullfile(subfolders{i}, 'cceps_files.m'), 'file')
                    warning('Excluding cceps_files.m from CCEP_Revised.');
                    continue; % Skip adding this folder to avoid conflicts
                end
                addpath(subfolders{i}); % Add other subfolders
            end
        end
    else
        error('Unknown CCEP version: %s', version);
    end
end