function [Pmf,angles] = spvm_pred(x,p)
% [Pmf,angles] = spvm_pred(x,p)
%    Predicts the probability mass function (Pmf) over discrete angle bins
%    using the Softplus von Mises model.
%
% Inputs:
%    x - Vector of predictor values.
%    p - A 2-element vector of model parameters.
%
% Outputs:
%    Pmf    - Matrix of predicted probability mass functions. Each row
%             corresponds to a distrinct element of x, and each column
%             corresponds to an angle bin.
%    angles - Vector of discrete angle values (in radians) used for the
%             prediction.
%
% Example:
%    [pmf, angles] = spvm_pred(x, [1, 0.05]);
%
% Subfunctions:
%    r2k   - Converts the vector of Softplus outputs to vector of von Mises
%            concentration parameters (k).
%    vmpdf - Computes a matrix of von Mises probability density functions.
%
% See also: spth_nll
%
x = x(:);
a = p(1);
b = p(2);
angles = (0:5).*(pi/3);
xPred = log(1+exp(x-a)).*b;
tPred = tanh(xPred);
kPred = r2k(tPred);
pdf = vmpdf(angles,0,kPred);
Pmf = pdf./sum(pdf,2);
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
for ii = 1:numel(r)
    r_ii = abs(r(ii));
    if (r_ii < 0.85) && (r_ii >= 0.53)
        k_ii = -0.4 + (1.39*r_ii) + (0.43/(1-r_ii));

    elseif r_ii < 0.53
        k_ii = (2*r_ii) + (r_ii^3) + (5*(r_ii^5)/6);

    else
        k_ii = ((r_ii^3) - (4*(r_ii^2)) + (3*r_ii))^-1;
    end
    k(ii) = k_ii * sign(r(ii));
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
% The pdf is normalised over all discrete angles in theta.
%
pdf = exp(k.*cos(theta-mu));
pdf = pdf ./ sum(pdf,2);
return