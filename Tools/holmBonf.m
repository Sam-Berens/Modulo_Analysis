function [rejectNull,thresholds] = holmBonf(pValues,alpha)

if nargin < 2
    alpha = 0.05;
end

m = numel(pValues);
[~,idx] = sort(pValues);
thresholds = nan(m,1);
for ii = 1:m
    thresholds(ii) = alpha / (m-ii+1);
end
thresholds(idx) = thresholds;
rejectNull = pValues < thresholds;
return