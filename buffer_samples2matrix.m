function [summary_outputs] = buffer_samples2matrix(x, samples_summarized)
%utility function to compute the running mean, mean amplitude deviation (MAD), and integration of the signal when having a vector and want to focus on signal segments

summary_outputs.x_buffered = buffer(x(1:end-rem(length(x),samples_summarized)), samples_summarized); % the vector in x is summarized in a matrix X with dimensions NxM, where N is the number of epochs and M is the number of samples within each epoch
summary_outputs.mean = mean(summary_outputs.x_buffered);
summary_outputs.MAD = mean(abs(summary_outputs.x_buffered - summary_outputs.mean)); % MeanAmplitudeDeviation
summary_outputs.sum_activity = sum(summary_outputs.x_buffered);


