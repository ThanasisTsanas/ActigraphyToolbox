function [Xsummarized] = summarize_matrix(x, p, q, function_used)
% Utility function to process a matrix
% Indicative call: [Xsummarized] = summarize_matrix(x, 5, 1, @median)

% Copyright (c) Athanasios Tsanas, 2014

% Last modified: 19 June 2014
%                16 May 2016, replacing 'data' with 'x' for easier debugging

%% Get inputs in convenient format
if nargin<4 || isempty(function_used)
    function_used = @mean; % default function for summarizing data every p rows
end

if nargin<3 || isempty(q)
    q = 1;
end

if nargin<2 || isempty(p)
    p = 5;
end

if (isvector(x))
   x = x(:); 
end

[N,M] = size(x);

%% Main processing
k1 = floor(N/p); k2 = M/q; % <-- k1 and k2 must evaluate to integers
[J,I] = meshgrid(p*(0:k1*k2-1)+k1*p*(q-1)*floor((0:k1*k2-1)/k1), (1:p*q)+(k1-1)*p*floor((0:p*q-1)/p));
J = round(J); p = floor(p);    

Xsummarized = reshape(function_used(reshape(x(I+J),p*q,k1*k2)),k1,k2);
