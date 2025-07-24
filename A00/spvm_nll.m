function [nll] = spvm_nll(b,x,c,Y)
% spvm_nll.m
% Sam Berens (s.berens@sussex.ac.uk)
% 22/07/2025
%
% Syntax:  nll = spvm_nll(b, x, c, Y)
%
% Description:
%    Computes the negative log-likelihood (nll) for the Softplus von Mises
%    model using a hyperbolic arc-tangent link function. Required inputs
%    are a set of model parameters, predictor values, and observed
%    responses (field indices).
%
% Inputs:
%    b   - A 2-element vector of model parameters.
%    x   - Vector of predictor values.
%    c   - Vector of target (correct) responses (as field indices).
%    Y   - Matrix of observed responses (field indices) for each trial.
%
% Outputs:
%    nll - The computed negative log-likelihood value.
%
% Example:
%    nll = spvm_nll([1, 0.5], x, c, Y);
%
% See also: spvm_pred
%
Theta = wrapTo2Pi((Y-repmat(c,1,6)).*(pi/3));
n = numel(x);
[pmf,angles] = spvm_pred(x,b);
nll = nan(n,1);
for iPred = 1:n
    p = pmf(iPred,:);
    t = Theta(iPred,:);
    t = t(~isnan(t));
    y = arrayfun(@(tt)find(abs(tt-angles)<1e-6),t);
    ty = y(1:(end-1));
    nll(iPred) = sum(...
        [0,log(1-cumsum(p(ty)))] - ...
        log(p(y)) ...
    );
end
nll(~isfinite(nll)) = -log(eps());
nll = sum(nll);
return