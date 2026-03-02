function S = invSoftplusLink()
%INVSOFTPLUSLINK  Custom inverse-softplus link for fitglme / fitglm*.
%
% This defines the link  g(mu) = softplus^{-1}(mu) = log(exp(mu) - 1)
% with inverse           g^{-1}(eta) = softplus(eta) = log(1 + exp(eta)).
%
% Domain notes:
%   - g(mu) is only defined for mu > 0.
%   - The inverse always returns mu > 0.
%
% Usage:
%   S = invSoftplusLink();
%   glme = fitglme(T, formula, 'Distribution', distName, 'Link', S);

    S = struct();
    S.Link             = @linkFun;
    S.Derivative       = @dlinkFun;
    S.SecondDerivative = @ddlinkFun;   % you may omit this field if allowed
    S.Inverse          = @invLinkFun;

end

% ---------- link: eta = g(mu) ----------
function eta = linkFun(mu)
    % eta = log(exp(mu) - 1) = log(expm1(mu))
    % mu must be > 0
    eta = log(expm1(mu));
end

% ---------- derivative: g'(mu) ----------
function d = dlinkFun(mu)
    % g'(mu) = exp(mu)/(exp(mu)-1) = exp(mu)/expm1(mu) = 1/(1-exp(-mu))
    d = exp(mu) ./ expm1(mu);
end

% ---------- second derivative: g''(mu) ----------
function dd = ddlinkFun(mu)
    % g''(mu) = -exp(mu)/(exp(mu)-1)^2 = -exp(mu)/expm1(mu)^2
    e1 = expm1(mu);
    dd = -exp(mu) ./ (e1.^2);
end

% ---------- inverse: mu = g^{-1}(eta) ----------
function mu = invLinkFun(eta)
    % mu = log(1 + exp(eta)) = softplus(eta)
    % Numerically stable implementation:
    mu = max(eta,0) + log1p(exp(-abs(eta)));
end
