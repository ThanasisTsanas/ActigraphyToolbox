function [data_compact, acc_measures] = acceleration_summary_measures(A, original_fs, new_fs, seconds_summarized)
%
% General call: [data_compact, acc_measures] = acceleration_summary_measures(A);
%               [data_compact, acc_measures] = acceleration_summary_measures(A, original_fs, new_fs, seconds_summarized);
%
%% Function to compute the four main acceleration summary measures:         ENMONZ, MAD, AI, ROCAM
%
% This function is based on my paper published in Sensors (2022)
%
% In general, read a binary file with three dimensional acceleration data using any standard funtions as required depending on the format originally saved
% for convenience, make sure the data is saved in the variable A -- see example below where I load data to understand format: should have three columns, each column indicates the x,y,z, axis, respectively 
%
% Aim: Compute the acceleration summary measures
%
% Inputs:  A       -> N by 3 matrix, N = observations. The columns MUST be
%                     in the form where the three columns indicate the
%                     acceleration in the x-, y-, and z-axis, respectively.
%                     The values in A indicate the acceleration measures in
%                     gravitational units (g). Therefore, the maximum value
%                     possible is dependent on the possible dynamic range
%                     in the device (typically +-8g).
%__________________________________________________________________________
% optional inputs:  
%     original_fs  -> the raw three dimensional acceleration data sample
%                     rate (e.g. 100 Hz in the CAPTURE-24 => put '100')      [default: 100]
%       new_fs     -> the new resampling frequency, now recommending to use 
%                     10 Hz for resampling data (efficient and fast)         [default: 10]
% seconds_summarized -> number of seconds over which we summarize the data
%                       to obtain e.g. minute-wise summary (with 60 sec)     [default: 60]
% =========================================================================
% Output:  
%     data_compact -> acceleration summary measures for specified segment
%                     (e.g. 60 seconds to obtain minute-wise assessments,
%                      as used in my papers)
%     acc_measures -> The acceleration measures (raw)
% =========================================================================
%
% Part of the "TsanasBOX" Toolbox
%
% -----------------------------------------------------------------------
% Useful references:
% 
% 1) A. Tsanas: Investigating wrist-based acceleration summary measures 
%    across different sample rates towards 24-hour physical activity and 
%    sleep profile assessment, Sensors, Vol. 22(16):6152, 2022
% 2) A. Tsanas, E. Woodward, A. Ehlers: Objective characterization of 
%    activity, sleep, and circadian rhythm patterns using a wrist-worn 
%    sensor: insights into post-traumatic stress disorder, JMIR mHealth and
%    uHealth, Vol. 8(4), pp. e14306, 2020 
% -----------------------------------------------------------------------
%
% Modification history
% --------------------
%  17 Aug 2022: Tidied function with full functionality for estimating the
%               four acceleration summary measures (ENMONZ, MAD, AI, ROCAM) 
%               from the raw three dimensional accelerometry data
%
% (c) Athanasios Tsanas
%
% ********************************************************************
% If you use this program please cite:
%
% 1) A. Tsanas: Investigating wrist-based acceleration summary measures 
%    across different sample rates towards 24-hour physical activity and 
%    sleep profile assessment, Sensors, Vol. 22(16):6152, 2022
% 2) A. Tsanas, E. Woodward, A. Ehlers: Objective characterization of 
%    activity, sleep, and circadian rhythm patterns using a wrist-worn 
%    sensor: insights into post-traumatic stress disorder, JMIR mHealth and
%    uHealth, Vol. 8(4), pp. e14306, 2020 
% ********************************************************************
%
% For any question, to report bugs, or just to say this was useful, email
% tsanasthanasis@gmail.com

%% Do some initial checks

% we MUST have an input matrix of size Nx3; generate an error otherwise
if(size(A,2)~=3)
    error('You need to provide a matrix which is in the form Nx3, i.e. contain three columns which correspond to the x,y,z axes!');
end

% check if the values presented in the matrix are not realistic, i.e.
% beyond 10g. For now just have it as a warning
if(max(abs(A(:)))>10)
    warning('Likely error, this indicates the acceleration was over 10g: check your dataset again!')
end

%% Set starting defaults

if nargin<4  || isempty(seconds_summarized)
    seconds_summarized = 60; %60; % summarize data every 60 seconds         RR!! have minute-level outputs
end

if nargin<3  || isempty(new_fs)
    new_fs = 10; % sampling frequency (Hz) for the processing (i.e. resample at this frequency == 10 Hz; could be 25, 50 or other values)
end

if nargin<2  || isempty(original_fs)
    original_fs = 100; % 100 Hz, this is the sample rate in CAPTURE-24
end

%% Resample accelerometry data to [10, 25, 50] Hz prior to any processing
data.xyz = resample(A, new_fs, original_fs); % subsample using FIR filter with anti-aliasing, also compensates for the delay introduced by the filter  

%% Now extract the different acceleration summary measures

acc_measures.EN = sqrt(data.xyz(:,1).^2 + data.xyz(:,2).^2 + data.xyz(:,3).^2); % this is the Gaussian computation, theoretically should be equal to gravity when not wearing watch/idle
acc_measures.ENMO = acc_measures.EN - 1;
acc_measures.ENMONZ = acc_measures.ENMO; acc_measures.ENMONZ(acc_measures.ENMONZ<0)=0; % ENMONZ (ENMO and making sure we have non-negative entries)
acc_measures.movement = [0;sqrt((diff(data.xyz(:,1))).^2 + (diff(data.xyz(:,2))).^2 + (diff(data.xyz(:,3))).^2)];

data_compact.ENMONZ = summarize_matrix(acc_measures.ENMONZ, new_fs*seconds_summarized, 1, @mean); 
[summary_outputs] = buffer_samples2matrix(acc_measures.EN, new_fs*seconds_summarized); data_compact.MAD = summary_outputs.MAD(:);

% the analysis for AI in the first place aims to find AI per second; then summing up to obtain the AI per minute
length_vector = length(acc_measures.EN); L = floor(length_vector/(new_fs));
x_buffered = buffer(data.xyz(1:end-rem(length_vector,L),1), new_fs); 
y_buffered = buffer(data.xyz(1:end-rem(length_vector,L),2), new_fs);
z_buffered = buffer(data.xyz(1:end-rem(length_vector,L),3), new_fs);

summary_buffered_mean_x = mean(x_buffered);
summary_buffered_mean_y = mean(y_buffered);
summary_buffered_mean_z = mean(z_buffered);

sigma_x = mean((x_buffered - summary_buffered_mean_x).^2);
sigma_y = mean((y_buffered - summary_buffered_mean_y).^2);
sigma_z = mean((z_buffered - summary_buffered_mean_z).^2);
AI_per_second = sqrt(max(mean([sigma_x;sigma_y;sigma_z],1),0));
remove_last_samples = rem(length(AI_per_second),seconds_summarized);
if(remove_last_samples>0), AI_per_second(length(AI_per_second)-remove_last_samples+1:end) = []; end
data_compact.AI = sum(buffer(AI_per_second, seconds_summarized))';

%ROCAM
% data_compact.ROCAM = moving_summary(summarize_matrix(acc_measures.movement, new_fs*seconds_summarized, 1, @mean), new_fs);
x = summarize_matrix(acc_measures.movement, new_fs*seconds_summarized, 1, @mean);
L = 10; data_compact.ROCAM = nanmean(buffer([NaN(1, floor(L/2)), x(:)', NaN(1, floor(L/2))], L, L-1, 'nodelay'))'; data_compact.ROCAM(length(x)+1:end) = [];

