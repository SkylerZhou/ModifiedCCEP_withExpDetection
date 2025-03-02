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

    % Add only the selected CCEP version
    if isfield(versions, version)
        addpath(genpath(versions.(version)));
        disp(['Using CCEP version: ', version]);
    else
        error('Unknown CCEP version: %s', version);
    end
end

