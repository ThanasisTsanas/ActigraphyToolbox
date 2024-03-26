function [PA_levels] = compute_PA_levels(data_compact)
%% Function to compute the Physical Activity (PA) levels from actigraphy data
%
% Indicative call: [PA_levels] = compute_PA_levels(data_compact);
%           ***********************************************
%
%% Compute various metrics from the Geneactiv data
%
% Aim: objectively characterize the different physical activity types
%
% Inputs:   data_compact     -> The Geneactiv data provided in compact form
%__________________________________________________________________________
% optional inputs:  
%   
% =========================================================================
% Outputs:  
%           PA_levels        -> structure with all the outputs for the
%                               different physical activity (PA) levels
%
% =========================================================================
%
% Part of the "Geneactiv processing" Toolbox by A. Tsanas
%
% -----------------------------------------------------------------------
% Useful references:
% 
%  [1]  A. Tsanas: Investigating wrist-based acceleration summary measures 
%       across different sample rates towards 24-hour physical activity and 
%       sleep profile assessment, Sensors, Vol. 22(16):6152, 2022

%  [2]  A. Tsanas, E. Woodward, A. Ehlers: Objective characterization of 
%       activity, sleep, and circadian rhythm patterns using a wrist-worn 
%       sensor: insights into post-traumatic stress disorder, JMIR mHealth 
%       and uHealth, Vol. 8(4), e14306, 2020 
%
%  [3]  V.T. van Hees, S. Sabia, K.N. Anderson, S.J. Denton, J. Oliver, 
%       M. Catt, J.G. Abell, M. Kivimäki, M.I. Trenell, A. Singh-Manoux: A 
%       novel, open access method to assess sleep duration using a wrist-worn 
%       accelerometer, PLoS One 10, e0142533, 2015 
%       (doi:10.1371/journal.pone.0142533)
%
% -----------------------------------------------------------------------
%
% Modification history
% ---------------------
% 22 Mar 2022: Function creation based on my Sensors 2022 paper
%
% Copyright (c) Athanasios Tsanas, 2020
%
% ********************************************************************
% If you use this program please cite:
%
%  [1] A. Tsanas: Investigating wrist-based acceleration summary measures 
%      across different sample rates towards 24-hour physical activity and 
%      sleep profile assessment, Sensors, Vol. 22(16):6152, 2022
%  [2] A. Tsanas, E. Woodward, A. Ehlers: Objective characterization of 
%      activity, sleep, and circadian rhythm patterns using a wrist-worn 
%      sensor: insights into post-traumatic stress disorder, JMIR mHealth 
%      and uHealth, Vol. 8(4), pp. e14306, 2020
%
% ********************************************************************
%
% For any question, to report bugs, or just to say this was useful, email
% tsanasthanasis@gmail.com

%% Check inputs and do computations accordingly for initializations

data_compact.xyz_movement = data_compact.ROCAM;
    
data_compact.dates = unique(floor(data_compact.time));
if(length(data_compact.dates)>1)
    remove_samples_begin = find(floor(data_compact.time) == floor(data_compact.time(1)), 1, 'last');
    remove_samples_end = find(floor(data_compact.time) == floor(data_compact.time(end)), 1, 'first');
else
    remove_samples_begin = [];
    remove_samples_end = [];
end

%remove prior and later samples (1st day and last day); then buffer the signal
xyz_movement = data_compact.xyz_movement; xyz_times = data_compact.time;
xyz_movement(remove_samples_end:end) = [];  xyz_times(remove_samples_end:end) = []; 
xyz_movement(1:remove_samples_begin) = []; xyz_times(1:remove_samples_begin) = [];
data_buffered = buffer(xyz_movement, 2*1440, 1440, 'nodelay');

% Create a matrix with the data only for a single day (i.e. using only
% the 1440 minutes on day, unlike the actogram-type of approach used
% for the matrix 'data_buffered' which has 2880 minutes data
for t = 1:size(data_buffered,2)
    data_buffered_diurnal(:,t) = data_buffered(1:1440,t); 
end
data_buffered_diurnal(:,t+1) = data_buffered(1441:end,t);

%% MAIN part of the function - compute the PA levels

PA_levels.dates = data_compact.dates(2:end-1);

% Follow the thresholds from the Tsanas (2022), Sensors paper. Note these are valid only when the sampling frequency was 10 Hz; will need to adapt according to the paper for other sample rates
a1 = (data_buffered_diurnal<0.035); % find sedentary activity
a2 = (data_buffered_diurnal<0.28 & data_buffered_diurnal>=0.035); % find light activity
a3 = (data_buffered_diurnal<0.6  & data_buffered_diurnal>=0.28); % find moderate activity
a4 = (data_buffered_diurnal>=0.6); % find vigorous activity
PA_levels.sedentary_activity = sum(a1); % sedentary activity: <1.5 METs
PA_levels.light_activity = sum(a2); % light activity:(1.5-4) METs
PA_levels.moderate_activity = sum(a3); % moderate activity: (4-7) METs
PA_levels.vigorous_activity = sum(a4); % vigorous activity: 7+ METs
PA_levels.sedentary_activity = PA_levels.sedentary_activity - NWTinfo.minutes_non_wear;

PA_levels.activity_map = a1 + 2.*a2 + 3.*a3 + 4.*a4;                         %RR!! 1 => sedentary, 2 => light, 3 => moderate, 4 => vigorous

if(plot_flag) % plot output
    PA_levels.figh = figure;
    bar(datetime(data_compact.dates(2:end-1), 'ConvertFrom', 'datenum'), [PA_levels.sedentary_activity; PA_levels.light_activity; PA_levels.moderate_activity; PA_levels.vigorous_activity], 'stacked');
    PA_levels.handle_PA_levels = gca;
    ylim([0 1440]); grid on; axis tight; PA_levels.handle_PA_levels.Children.XTickLabelRotation = 90;
    xlabel('Dates', 'FontSize', 16, 'FontWeight','bold'); ylabel('Minutes', 'FontSize', 16, 'FontWeight','bold'); 
    title(['Activity types: ', datestr(data_compact.dates(2), 'dd mmm yyyy'), ' - ', datestr(data_compact.dates(end), 'dd mmm yyyy')], 'FontWeight', 'Bold', 'FontSize', 14);
    lh = legend({'Sedentary', 'Light', 'Moderate', 'Vigorous'}, 'Location', 'southeast'); %lh.Orientation = 'horizontal'; lh.NumColumns = 5; lh.Location = 'south'; lh.FontSize = 8; 
    newcolors = [colormap(parula(4)); 1 0.6 0.4]; 
    colororder(newcolors)
end
