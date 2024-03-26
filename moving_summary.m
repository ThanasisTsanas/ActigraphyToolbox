function y = moving_summary(x, L, summarizing_operator, overlap_samples, operation_flag)
% Function for summarizing a vector based on the summarizing operator

% (c) Copyright, A. Tsanas 2015

% Modification history
% ---------------------
%  5 Aug 2015: function creation
% 22 Jul 2017: introduced 'operation_flag'

%% Get initial defaults
if(nargin<2 ||isempty(L))
    L = 5; % by default apply to five successive samples
end

if(nargin<3 ||isempty(summarizing_operator))
    summarizing_operator = @nanmean;
end

if(nargin<4 ||isempty(overlap_samples))
    overlap_samples = L-1;
end

if(nargin<5 ||isempty(overlap_samples))
    operation_flag = 'centered';
end

%% Main part of the function
switch(operation_flag)
    case {'centered', 1}
        X = buffer([NaN(1, floor(L/2)), x(:)', NaN(1, floor(L/2))], L, overlap_samples, 'nodelay');
    case {'start', 'first', 2}
        X = buffer([x(:)', NaN(1, floor(L))], L, overlap_samples, 'nodelay');
    otherwise
        disp('Not correctly set mode operation for the mask!');
end
y = summarizing_operator(X); 

if(overlap_samples~=0)
    y(length(x)+1:end) = [];
end

y=y';
