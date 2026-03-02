function [lq] = BN_ll(b,x,theta)
%takes b, a matrix of model paremeter paremeters dims: nSamples x 2
%takes x, a predictor value
%takes resp, a
%outputs vector of loglikelihoods for all samples

[pmf,bin] = BN_pred(x,b);

theta = theta(~isnan(theta)); %chopp off the nan attempts
nTry = size(theta,1);
norm = 1;
lq = 0;
for k = 1:nTry
    t = theta(k);
    isCorr = abs(t)<1e-6;
    [~,idx] = ismember(isCorr,bin);
    linidx = sub2ind(size(pmf),(1:size(pmf,1))',idx);
    cP = pmf(linidx); %current likelihood of model given the observation / probability of data given the model
    prInco = pmf(:,1);
    if k>1
        norm = 1 - ((k-1).*prInco/5);
    end
    lq = lq + log(cP) - log(norm);  %sum across attmepts to get log of pr/likelihood for each trial
end


if ~isfinite(lq)
    warning('Non-finite log-likelihood detected');
end

return