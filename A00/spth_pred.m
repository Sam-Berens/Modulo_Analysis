function [pmf,angles] = spth_pred(x,b)
% spth_pred.m
% Sam Berens (s.berens@sussex.ac.uk)
% 24/03/2025
%
% Syntax:  [pmf, angles] = spth_pred(x, b)
%
% Description:
%    Predicts the probability mass function (pmf) over discrete angle bins
%    using a Softplus-tanh model with the von Mises distribution. The
%    function computes an intermediate prediction, converts it to a
%    concentration parameter using r2k, and then calculates a von Mises
%    probability density function that is normalized to yield the pmf.
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
%    [pmf, angles] = spth_pred(x, [1, 0.5]);
%
% Subfunctions:
%    r2k   - Converts the Softplus-tanh output to a von Mises concentration
%            parameter k.
%    vmpdf - Computes the von Mises probability density function.
%
% See also: spth_nll
%

angles = (0:5).*(pi/3);
xPred = log(1+exp(x-b(1)))*b(2);
tPred = tanh(xPred);
kPred = r2k(tPred);
pdf = vmpdf(angles,0,kPred);
pmf = pdf./sum(pdf,2);
return

function [k] = r2k(r)
% r2k converts the Softplus-tanh output values to a von Mises concentration
% parameter k.
%
% Inputs:
%    r - A vector of Softplus-tanh transformed predictor values.
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
% vmpdf computes the von Mises probability density function.
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
% The pdf is normalized by the factor c, which involves the zeroth-order
% Bessel function.
%

c = 1./(2.*pi.*besseli(0,k));
pdf = c .* exp(k.*cos(theta-mu));
return