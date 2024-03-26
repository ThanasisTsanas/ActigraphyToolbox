function [diurnal_measures] = compute_diurnal_measures(data_compact)
%% Function to compute various measures from the actigraphy data
%
% Indicative call: [diurnal_measures] = compute_diurnal_measures(data_compact);
%           ***********************************************
%
% Aim: objectively characterize the Geneactiv data extracting diurnal measures
%
% Inputs:   data_compact     -> The Geneactiv data provided in compact form
%__________________________________________________________________________
% optional inputs:  
%   
%               N/A
% =========================================================================
% Outputs:  
%           diurnal_measures  -> structure with all the outputs. Each struct
%                               entry is a vector which is the same as the
%                               number of days in the data with the
%                               exception of the first and last days
%                               (because these would be incomplete)
%
% =========================================================================
%
% Part of the "Actigraphy Analysis Toolbox" by A. Tsanas
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
%
% -----------------------------------------------------------------------
%
% Modification history
% ---------------------
% 18 Mar 2024: code cleaning, preparing the function for open-sourcing
%
% Copyright (c) Athanasios Tsanas, 2024
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

%% Check inputs and do computations accordingly for initializations

data_compact.xyz_movement = data_compact.ROCAM;
if(length(data_compact.ROCAM)>1)
    remove_samples_begin = find(floor(data_compact.time) == floor(data_compact.time(1)), 1, 'last');
    remove_samples_end = find(floor(data_compact.time) == floor(data_compact.time(end)), 1, 'first');
else
    remove_samples_begin = [];
    remove_samples_end = [];
end
xyz_movement = data_compact.xyz_movement; xyz_times = data_compact.time;
xyz_movement(remove_samples_end:end) = [];  xyz_times(remove_samples_end:end) = []; 
xyz_movement(1:remove_samples_begin) = []; xyz_times(1:remove_samples_begin) = [];
data_buffered = buffer(xyz_movement, 2*1440, 1440, 'nodelay');
for t = 1:size(data_buffered,2)
    data_buffered_diurnal(:,t) = data_buffered(1:1440,t); 
end
data_buffered_diurnal(:,t+1) = data_buffered(1441:end,t);

%% Main part
diurnal_measures.dates = unique(floor(xyz_times));
mean_overall_activity = mean(nanmean(data_buffered_diurnal));
[Xsummarized] = summarize_matrix(data_buffered_diurnal, 60, 1, @mean); 
for t = 1:size(data_buffered_diurnal,2)
    y_M10(:,t) = moving_summary(data_buffered_diurnal(:,t), 10*60, [], [], 'start');
    y_L5(:,t) = moving_summary(data_buffered_diurnal(:,t), 5*60, [], [], 'start');
    Xsummarized_30min_overlap(:,t) = moving_summary(data_buffered_diurnal(:,t), 1*60, [], 1*30, 'start');
end
[diurnal_measures.M10val, diurnal_measures.M10idx] = max(y_M10); 
[diurnal_measures.L5val, diurnal_measures.L5idx] = min(y_L5(1:1000,:));

diurnal_measures.M10idx(diurnal_measures.M10idx<0) = 1 - diurnal_measures.M10idx(diurnal_measures.M10idx<0);
diurnal_measures.L5idx(diurnal_measures.M10idx<0) = 1 - diurnal_measures.L5idx(diurnal_measures.L5idx<0);
diurnal_measures.RAval = (diurnal_measures.M10val - diurnal_measures.L5val)./(diurnal_measures.M10val + diurnal_measures.L5val);

Xsummarized_sorted = sort(Xsummarized);
diurnal_measures.IS1 = (mean((Xsummarized - mean_overall_activity).^2))./(mean((data_buffered_diurnal - mean_overall_activity).^2));
diurnal_measures.IS2 = (nanmean((Xsummarized_30min_overlap - mean_overall_activity).^2))./(nanmean((data_buffered_diurnal - mean_overall_activity).^2));
diurnal_measures.IV1 = mean(diff(Xsummarized).^2)./mean((Xsummarized - repmat(mean(Xsummarized), size(Xsummarized,1), 1)).^2); 
diurnal_measures.IV2 = mean(diff(data_buffered_diurnal).^2)./mean((data_buffered_diurnal - repmat(mean(data_buffered_diurnal), size(data_buffered_diurnal,1), 1)).^2); %1440 minutes
diurnal_measures.IV3 = nanmean(diff(Xsummarized_30min_overlap).^2)./(nanmean((Xsummarized_30min_overlap - mean_overall_activity).^2)); 
diurnal_measures.ratioIV = diurnal_measures.IV1./diurnal_measures.IV2;
diurnal_measures.L5 = mean(Xsummarized_sorted(1:5,:));
diurnal_measures.M10 = mean(Xsummarized_sorted(end-9:end,:));
diurnal_measures.RA = (diurnal_measures.M10 - diurnal_measures.L5)./(diurnal_measures.M10 + diurnal_measures.L5);
diurnal_measures.activity_percentiles5_25_50_75_95 = prctile(Xsummarized, [5 25 50 75 95]);
diurnal_measures.activity_mean = prctile(data_buffered_diurnal, [5 25 50 75 95]);
diurnal_measures.activity_TKEO = mean(TKEO(data_buffered_diurnal));
diurnal_measures.activity_RMSSD = RMSSD(data_buffered_diurnal);
diurnal_measures.mean_activity = mean(data_buffered_diurnal);

[diurnal_measures.PA_levels] = compute_PA_levels(data_compact);
[diurnal_measures.multiscale_entropy_data] = compute_MSE(data_compact); 

%% Sleep-based characteristics
% will need to clean up my sleep code and open-source here
% as a proxy for now, can use the ROCAM-sleep thresholded outputs



end %% End of main function

%% =====================================================================

function [multiscale_entropy_data] = compute_MSE(data_compact)   
    variable_diurnal = 'ROCAM';
    data_compact.xyz_movement = moving_summary(mean(data_compact.(variable_diurnal),2), 5, @nanmean);

    if(length(data_compact.dates)>1)
        remove_samples_begin = find(floor(data_compact.time) == floor(data_compact.time(1)), 1, 'last');
        remove_samples_end = find(floor(data_compact.time) == floor(data_compact.time(end)), 1, 'first');
    else
        remove_samples_begin = [];
        remove_samples_end = [];
    end
    
    xyz_movement = data_compact.xyz_movement; xyz_times = data_compact.time;
    xyz_movement(remove_samples_end:end) = [];  xyz_times(remove_samples_end:end) = []; 
    xyz_movement(1:remove_samples_begin) = []; xyz_times(1:remove_samples_begin) = [];
    data_buffered = buffer(xyz_movement, 2*1440, 1440, 'nodelay');
    
    for t = 1:size(data_buffered,2)
        data_buffered_diurnal(:,t) = data_buffered(1:1440,t); 
    end
    data_buffered_diurnal(:,t+1) = data_buffered(1441:end,t);
    sanity_check = cellstr(datestr(xyz_times)); % just for me to ensure we have days starting and finishing at 00:00
    composite_multiscale_entropy_data_all = CMSE(xyz_movement,1440);

    for t = 1:size(data_buffered_diurnal,2)
        clear composite_multiscale_entropy_data multiscale_entropy_data1
        composite_multiscale_entropy_data = CMSE(data_buffered_diurnal(:,t), 120);
        CMSE_results(:,t) = composite_multiscale_entropy_data([5 10 30 60]);
    end

    multiscale_entropy_data.CMSE = [repmat(composite_multiscale_entropy_data_all([5 10 30 60 120]), 1, size(data_buffered_diurnal,2)); CMSE_results];
end % end of compute_MSE
