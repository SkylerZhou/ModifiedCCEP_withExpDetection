function locations = cceps_files_example

% dir related to ieeg
locations.pwfile = "dir to your ieeg password file";
locations.ieeg_folder = 'dir to your IEEGToolbox';
locations.loginname = 'type your ieeg login name here';

% dir to OUTPUT file from 1st versions
locations.results_folder = 'dir to store the output of the 1st pipeline (erin version). i have my first-version results stored at different directories so i have both result_folder and firstOut_dir, if you have them stored at the same place these two would have the same dir.'; 
locations.firstOut_dir = 'dir to store the output of the 1st pipeline (erin version)';

% dir to INPUT scripts & pt data & eletrode dist info for 3rd version
locations.script_folder = 'dir to your ccep detection script';
locations.data_folder = 'dir to the pt_mat file';
locations.mount_dir = 'dir to your electrode coordinate files';

% dir to OUTPUT for 3rd version
locations.thirdOut_dir = 'dir to store the output of the 3rd pipeline (skyler version)';

end