function locations = cceps_files

% dir related to ieeg
locations.pwfile = "/Users/zhouzican/Documents/MATLAB/master_thesis/ieeg_login/sky_ieeglogin.bin";
locations.ieeg_folder = '/Users/zhouzican/Documents/MATLAB/master_thesis/ieeg-matlab-1.14.60/IEEGToolbox/';
locations.loginname = 'skylerz';

% dir to OUTPUT file from 1st versions
locations.results_folder = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP/ccep_result/'; 
% patientOriOut_dir, should modify to the dir under ccep_revised later
locations.firstOut_dir = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP_Revised/ccep_result/first_pipeline/';
% locations.firstOut_dir = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP/ccep_result/ori_pipeline';

% dir to INPUT scripts & patient info & eletrode dist data for 3rd version
locations.script_folder = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP_ExpDetection/ccep_script/';
locations.data_folder = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP_ExpDetection/pt_mat/';
locations.mount_dir = '/Volumes/USERS/skylerz/coor_data/coor_files/'; % elec distance data

% dir to OUTPUT for 3rd version
%locations.patientNewOut = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP/ccep_result/new_pipeline_keptonly';
% patientNewOut
locations.thirdOut_dir = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP_ExpDetection/ccep_result/third_pipeline_threshAmp4.5/';
% locations.thirdOut_dir = '/Users/zhouzican/Documents/MATLAB/master_thesis/CCEP_ExpDetection/ccep_result/third_pipeline';

% files containing electrodes distance info
%locations.coord_folder_dir1 = '/mnt/leif/littlab/data/Human_Data/CNT_iEEG_BIDS/';
%locations.coord_folder_dir2 = '/mnt/leif/littlab/data/Human_Data/recon/BIDS_penn/';
%locations.coord_file_dir = '/derivatives/ieeg_recon/module3/';

% dir for hypothesis testing
locations.sz_ccep_overlap = '/Users/zhouzican/Documents/CNT/DSOSD/SZ_hypo1_data/goldStandard/sz_ccep_overlap.csv';

end