%% This file aims to plot EC's validations using RW's format 

% load ori_out 
ptT = readtable(['/Users/zhouzican/Documents/MATLAB/toolboxs/CCEP/pt_mat/','master_pt_list.xlsx']);
patient_files = string(strcat(ptT.HUPID, '.mat'));    
start_patient = 40;
num_patient = 40;


for n = start_patient:num_patient
    patient_file = fullfile('toolboxs', 'CCEP', 'ccep_result', 'new_pipeline', patient_files(n));
    out = load(patient_file);
    ori_out = out.pt_out;

    
end
