function [pmf,angles] = spvm_pred(x,b)
% spvm_pred.m
% Sam Berens (s.berens@sussex.ac.uk)
% 24/03/2025
%
% Syntax:  [pmf, angles] = spvm_pred(x, b)
%
% Description:
%    Predicts the probability mass function (pmf) over discrete angle bins
%    using the Softplus von Mises model with a hyperbolic arc-tangent link
%    function. The function computes an intermediate prediction, converts
%    it to a concentration parameter using r2k, and then calculates a von
%    Mises probability density function that is normalised to yield the
%    pmf.
%
% Inputs:
%    x - Vector of predictor values.
%    b - A 2-element vector of model parameters.
%
% Outputs:
%    pmf    - Matrix of predicted probability mass functions for each 
%             trial. Each row corresponds to a trial, and each column 
%             corresponds to an angle bin.
%    angles - Vector of discrete angle values (in radians) used for the 
%             prediction.
%
% Example:
%    [pmf, angles] = spvm_pred(x, [1, 0.5]);
%
% Subfunctions:
%    r2k   - Converts the inverse link function output to a von Mises
%            concentration parameter, kappa.
%    vmpdf - Computes the von Mises probability density function.
%
% See also: spvm_nll
%
if numel(b) ~= 2
    error('b must be a 2-vector');
end
angles = (0:5).*(pi/3);
xPred = log(1+exp(x-b(1)))*b(2);
tPred = tanh(xPred);
kPred = r2k(tPred);
pdf = vmpdf(angles,0,kPred);
pmf = pdf./sum(pdf,2);
return

function [k] = r2k(r)
% r2k converts the inverse link function output values to a von Mises
% concentration parameter, kappa.
%
% Inputs:
%    r - A vector of Softplus tanh transformed predictor values.
%
% Outputs:
%    k - A vector of concentration parameters corresponding to r.
%
% The conversion uses piecewise definitions based on abs(r).
%

k = nan(size(r));
for ii = 1:1:numel(r)
    cr = abs(r(ii));
    if (cr < 0.85) && (cr >= 0.53)
        ck = -0.4 + (1.39*cr) + (0.43/(1-cr));
    elseif cr < 0.53
        ck = (2*cr) + (cr^3) + (5*(cr^5)/6);
    else
        ck = ((cr^3) - (4*(cr^2)) + (3*cr))^-1;
    end
    k(ii) = ck*sign(r(ii));
end
k(k>700) = 700;
k(k<-700) = -700;
return

function [pdf] = vmpdf(theta,mu,k)
% vmpdf computes the (unnormalised) von Mises PDF.
% ... unnormalised, as it is normed into a valid PMF above.
%
% Inputs:
%    theta - Vector of angle values (in radians).
%    mu    - Mean direction (in radians).
%    k     - Concentration parameter.
%
% Outputs:
%    pdf   - Vector of probability density values computed at each angle in
%            theta.
%

% c = 1./(2.*pi.*besseli(0,k));
pdf = exp(k.*cos(theta-mu));
return
