function [lq] = spvm_ll(P,x,theta)
%MODIFIED VERSION of spth_nll.m
% Description:
%    Computes the negative log-likelihood (nll) for a Softplus-tanh model
%    using the von Mises distribution given a set of model parameters,
%    predictor values, and observed angle responses.
%
% Inputs:
%    b     - A 2-element vector of model parameters.
%    x     - Vector of predictor values.
%    theta - Matrix of observed angle responses (in radians) for each 
%            trial.
% Outputs:
%    lq   - Vector of log-likelihoods for all samples 

[pmf,angles] = spvm_pred(P,x);
theta = theta(~isnan(theta));
nTry = numel(theta);
massExcluded = 0;
lq = 0;
for k = 1:nTry
    t = theta(k);
    sResp = abs(t-squeeze(angles)') < 1e-6;
    p = pmf(:,sResp);
    lq = lq + log(p) - log(1-massExcluded);
    massExcluded = massExcluded + p;
end
if ~isfinite(lq)
    warning('Non-finite log-likelihood detected');
end
return