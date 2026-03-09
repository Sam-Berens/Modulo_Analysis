function [Pmf,bin] = spbn_pred(x,p)
% [Pmf,bin] = spbn_pred(x,b)
%    Predicts the probability mass function (Pmf) using the Softplus
%    binomial model.
x = x(:);
a = p(1);
b = p(2);
xPred = log(1 + exp(x-a)) .* b;
prCorr = 1./(1+exp(log(5)-xPred));
prInco = 1 - prCorr;
Pmf = [prCorr,prInco];
bin = [1,0];
return