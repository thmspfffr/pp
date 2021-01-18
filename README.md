# pp

# Prepare data

1) pp_prepare_data.m

# Analyse data

combines MEG, eye tracker recordings as well as MEG channel information into one dataset for both resting state and task

2a) pp_src_pupil_power_correlations.m
  - computes pupil-MEG correlations during rest, hamburg data set only
2b) pp_gla_src_pupil_power_correlations.m
  - computes pupil-MEG correlations during rest, glasgow data set only
2c) pp_mue_src_pupil_power_correlations.m
  - computes pupil-MEG correlations during rest, muenster data set only
  
3a) pp_hh_src_powerspectra.m
3b) pp_gla_src_powerspectra.m
3c) pp_mue_src_powerspectra.m
  - compute power spectra in segments: result used for fitting fooof model (see 4)
  - separately for HH, GLA and MUE
  
4a) pp_src_hh_fooof.py
  - fit fooof model to results obtained from (3)
4b) pp_src_gla_fooof.py
4c) pp_src_mue_fooof.py

5) pp_cnt_src_pupil_power_correlations.m
  - Compute pupil-MEG correlations during task (hamburg only)
6) pp_hh_task_src_powerspectra.m
7) pp_src_hh_task_fooof.m

# Plotting

figure1.m
figure2.m

etc. (to be updated - 14/01/2020)




  
