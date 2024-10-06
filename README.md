# ModifiedCCEP_withExpDetection

## 1. to_run_newPipeline
Updated on three aspects based on https://github.com/erinconrad/CCEPS/tree/main and https://github.com/RudyWh/Modified_CCEP_Detection to fit new cohort:
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
