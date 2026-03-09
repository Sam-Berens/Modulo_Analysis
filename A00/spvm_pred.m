function [Pmf,angles] = spvm_pred(P,x)
% [Pmf,angles] = spvm_pred(P,x)
%    Predicts the probability mass function (Pmf) over discrete angle bins
%    using the Softplus von Mises model.
%
% Inputs:
%    x - Vector of predictor values.
%    p - A 2-element vector of model parameters.
%
% Outputs:
%    Pmf    - Matrix of predicted probability mass functions. Each row
%             corresponds to a disct row in P, each column corresponds to a
%             distrinct element of x, and each page corresponds to an angle
%             bin.
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
if size(P,2) ~= 2
    error('P must be an n-by-2 matrix');
end
x = x(:)';
angles = reshape((0:5).*(pi/3),1,1,[]);

% nSamples = size(P,1);
% nTrials = numel(x);
% Pmf = nan(nSamples,nTrials,6);

a = P(:,1);
b = P(:,2);
xPred = log(1+exp(x-a)).*b;
tPred = tanh(xPred);
kPred = r2k(tPred);
Pmf = vmPmf(angles,0,kPred);
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

function [Pmf] = vmPmf(theta,mu,k)
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
% The pmf is normalised over all discrete angles in theta.
%
Pdf_ish = exp(k.*cos(theta-mu));
Pmf = Pdf_ish ./ sum(Pdf_ish,3);
return