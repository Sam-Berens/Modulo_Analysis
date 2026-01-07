function [lq] = vM_ll(b,x,theta)
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

[pmf,angles] = vM_pred(x,b);
theta = theta(~isnan(theta)); %chop off the nan attempts
nTry = size(theta,1);
cumP = 0; % cumulative probability 
lq = 0;
for k = 1:nTry
    t = theta(k);
    sResp = abs(t-angles)<1e-6;
    cP = pmf(:,sResp); %current probability of response given the model
    if k==1
        lq = log(cP);
    else
        condP = cP ./ (1 - cumP); %conditional probability for each attempt
        lq = lq + log(condP);
    end
    cumP = cumP + cP;

end
lq = sum(lq,2); %dim 2 is attempts dim 1 is samps
if ~isfinite(lq)
    warning('Non-finite log-likelihood detected');
end

return
