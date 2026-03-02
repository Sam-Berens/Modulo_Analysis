function [alphas,rejectNull] = holmBonf(alpha, ps)

m = numel(ps);
[ps,idx] = sort(ps);
alphas = nan(m,1);
rejectNull = nan(m,1);
for ii=1:m
    alphas(ii) = alpha / (m-ii+1);
    rejectNull(ii) = ps(ii) < alphas(ii);
end

alphas(idx) = alphas;
rejectNull(idx) = rejectNull;

return