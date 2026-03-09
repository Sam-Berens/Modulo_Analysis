function [lq] = spbn_ll(P,x,theta)
%takes b, a matrix of model paremeter paremeters dims: nSamples x 2
%takes x, a predictor value
%takes resp, a
%outputs vector of loglikelihoods for all samples










[pmf,bin] = spbn_pred(P,x);
theta = theta(~isnan(theta)); 
nTry = size(theta,1);
massExcluded = 0;
lq = 0;
for k = 1:nTry
    t = theta(k);
    isCorrect = double(abs(t)<1e-6);
    [~,idxChosen] = ismember(isCorrect,bin);
    p = pmf(:,idxChosen);
    lq = lq + log(p) - log(1-massExcluded);
    massExcluded = massExcluded + p;
end
if ~isfinite(lq)
    warning('Non-finite log-likelihood detected');
end
return