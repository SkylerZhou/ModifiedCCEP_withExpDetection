# ModifiedCCEP_withExpDetection
Note: Opening multiple CCEP folders might result in conflict of functions. Remember to remove one folder from path if that is the case.

## 0. run setup_ccep.m (optional)
If there is frequent need to switch between different versions of the CCEP detectors, can run this script using setup_ccep('specific_version') to avoid cross-using of the directories in cceps_files.m
1. First version: https://github.com/erinconrad/CCEPS
2. Second version: https://github.com/RudyWh/Modified_CCEP_Detectio
3. Thrid version: https://github.com/SkylerZhou/ModifiedCCEP_withExpDetection

## 1. to_run_newPipeline
Updated on three aspects based on https://github.com/erinconrad/CCEPS/tree/main and https://github.com/RudyWh/Modified_CCEP_Detection to fit new cohort.
Modifications include:
1. Added an **exp** component in AA_new_build_network.m, which utilize Curve Fitting Toolbox to fit an *standardized & low-pass filtered* exponential function to approximately 0-0.3 second of the recording. Reject the signal if its goodness-of-fit is bigger than 0.6. 
2. Lower the **thresh_amp** from 0.65 to 0.6.
3. Delete the **deriv** component in AA_alternative_filtering.m and AA_new_build_network.m.

## 2. validation_on_newPipeline
Supplementary files for graphing and checking whether expoential waves and stim artifacts are correctly rejected. 
1. SZ_validateExp.m randomly select 25 signals that were rejected due to the exp component. The signals were mapped out so that we could visually examine if they were actually of expoential shape and thus correctly rejected. 
2. Stim artifacts noticeably increase among the 25 randomly selected kept figures after processing using the new pipleine. To aid in checking if these artifacts were retained in both original and new pipeline, SZ_validateStimArtifact.m were built to compare them in both versions of the pipeline using the same set of coordinates. 

## 3. sanity_checks
Conducted pair-wise comparison between N1&N2 amplitudes and latencies, as well as electrode distances as a simple sanity check on the new pipeline. 
1. SZ_adjust_network_to_remove_rejects.m convert amplitude and latency of the N1&N2 that were rejected to NaN.
2. SZ_add_elecs_distance.m adds a electrode distance matrix. 
3. SZ_spearman_ampLatDist_ptLevel.m compute pairwise correlation between amp, lat, and dist at the patient level for sanity check (we expect amp to be negatively correlated with latency and distance, while latency and distance to be postively correlated). 

## 4. convert2_csv
1. SZ_keptOnly.m add the eletrode eucleadian distance matrix to out.other.elecs_dist; creat out.elecs.n1_adj and out.elecs.n2_adj to store the n1&n2 amplitude, latency, and elec distance info of the signals that are retained after the CCEP detector. 
2. SZ_mat2csv.m convert the mat output to csv for downstream python processing 