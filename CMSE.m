function E = CMSE(data,scale)

%Composite Multiscale Entropy Algorithm

% See: "Time Series Analysis Using Composite Multiscale Entropy", Entropy 2013, 15, 1069-1084; doi:10.3390/e15031069

if nargin<2 || isempty(scale)
    scale = 20; % following the paper's recommendation
end

r = 0.15*std(data);

progressbar;
E = zeros(scale,1);
for i = 1:scale % i:scale index
    for j = 1:i % j:croasegrain series index
        buf = CoarseGrain(data(j:end),i);
        E(i) = E(i)+ SampEn(buf,r)/i;
    end
    progressbar(i/scale) % Update figure
%     i
end

%Coarse Grain Procedure. See Equation (2)
% iSig: input signal ; s : scale numbers ; oSig: output signal
function oSig=CoarseGrain(iSig,s)
N=length(iSig); %length of input signal
for i=1:1:N/s
    oSig(i)=mean(iSig((i-1)*s+1:i*s));
end

%function to calculate sample entropy. See Algorithm 1
function entropy = SampEn(data,r)
l = length(data);
Nn = 0;
Nd = 0;

for i = 1:l-2
    for j = i+1:l-2
        if abs(data(i)-data(j))<r && abs(data(i+1)-data(j+1))<r
            Nn = Nn+1;
            if abs(data(i+2)-data(j+2))<r
                Nd = Nd+1;
            end
        end
    end
end

entropy = -log(Nd/Nn);
