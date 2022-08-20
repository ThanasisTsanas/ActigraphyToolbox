# ActigraphyToolbox

This is a collection of MATLAB files to facilitate the analysis of actigraphy data.

****************************************
I intend to be uploading different functions as I publish in the area of actigraphy data processing.

****************************************
For now I have included some sample data (see sample_data.mat) and a demo to compute the different acceleration summary measures, as described in my Sensors 2022 paper. Specifically, I compute ENMONZ, MAD, AI, ROCAM (the latter being the novel acceleration summary measure I proposed in the paper). The outputs of these acceleration summary measures can be used to compute sleep and different levels of Physical Activity (PA) - see my Sensors 2022 paper for details.

ROCAM builds on the ideas I have developed previously in my JMIR mHealth 2020 paper on the use of 'movement' as an acceleration summary measure.

Run the file 'demo_acceleration_summary_measures.m'

****************************************

Copyright (c) Athanasios Tsanas, 2022

**If you use this program please cite:**

1) A. Tsanas: Investigating wrist-based acceleration summary measures 
   across different sample rates towards 24-hour physical activity and 
   sleep profile assessment, Sensors, Vol. 22(16):6152, 2022

2) A. Tsanas, E. Woodward, A. Ehlers: Objective characterization of 
  activity, sleep, and circadian rhythm patterns using a wrist-worn 
  sensor: insights into post-traumatic stress disorder, JMIR mHealth and
  uHealth, Vol. 8(4), pp. e14306, 2020 
