
% Load indicative three-dimensional acceleration data which is in the form Nx3
load('sample_data.mat') % once loaded, the data appears in the matrix A (which is a matrix of Nx3 dimensions)

original_fs = 100; % the CAPTURE-24 data has been sampled at 100 Hz, so use that here for this example
new_fs = 10; % as I have shown in my Sensors2022 paper, 10 Hz is broadly sufficient for sleep and Physical Activity (PA) assessment
seconds_summarized = 60; % get overall estimates every 60 seconds, i.e. have minute-by-minute overall estimate of the acceleration measures which will be matched to minute-wise actigrpahy labels (as we have in CAPTURE-24, see my paper in Sensors2022)

% call the function to compute the acceleration summary measures
[data_compact, acc_measures] = acceleration_summary_measures(A, original_fs, new_fs, seconds_summarized);

% the output of interest here is in the variable 'data_compact', which
% contains the four acceleration summary measures used in my paper



