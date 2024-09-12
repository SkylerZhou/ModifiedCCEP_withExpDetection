# ModifiedCCEP_withExpDetection

## Updated on three aspects of https://github.com/RudyWh/Modified_CCEP_Detection to fit our cohort:
1. Added an **exp** component in AA_new_build_network.m, which utilize Curve Fitting Toolbox to fit an exponential function to approximately 0-0.3 second of the recording. Reject the signal if its goodness-of-fit is bigger than 0.7.
2. Lower the **thresh_amp** from 0.65 to 0.6.
3. Delete the **deriv** component in AA_alternative_filtering.m and AA_new_build_network.m.
